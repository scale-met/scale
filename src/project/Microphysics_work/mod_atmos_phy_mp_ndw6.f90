!-------------------------------------------------------------------------------
!
!+  NICAM Double moment Water 6 scheme
!
!-------------------------------------------------------------------------------
module mod_atmos_phy_mp
  !-----------------------------------------------------------------------------
  !
  !++ Description: 
  !       This module contains subroutines for the ndw6 parametrization.
  !
  !       
  !++ Current Corresponding Author : T.Seiki
  ! 
  !++ History: NDW6
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
  !                   mod_mp_nsw6.f90 in NICAM
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_stdio, only: &
     IO_FID_LOG,  &
     IO_L
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
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: mp_ndw6_init
  private :: mp_ndw6
  private :: terminal_velocity     
  private :: mp_ndw6_diag_volume           
  private :: mp_ndw6_effective_radius      
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters
  !
  integer, parameter :: HYDRO_MAX  = 6    ! total number of mixing ratio of water 
  character(len=32), save :: WLABEL(11)

  ! empirical value from Meyers etal.(1991), 1[/liter] = 1.d3[/m3]
  real(8), private, parameter :: nqmin(6) = (/ 0.D0, 1.D4, 1.D0, 1.D0, 1.D-4, 1.D-4 /) ! [1/m3]
  ! refer to Seifert(2002) (Dr. Thesis, Table.5.1)
  ! max mass, for D_min=79um, 2mm, 5mm, 1cm, 1cm
  real(8), private, parameter :: xqmax(6) = (/ 0.D0, 2.6D-10, 5.0D-6, 1.377D-6, 7.519D-6, 4.90D-5 /)! [kg]
  ! SB06, Table 1.
  ! min mass, for D_min=2um, 79um, 10um, 20um, 100um
  real(8), private, parameter :: xqmin(6) = (/ 0.D0, 4.20D-15, 2.60D-10, 3.382D-13, 1.847D-12, 1.230D-10 /)! [kg]



  ! for all processes
  ! SB06, Table 1.
  real(8), parameter :: xc_min = 4.20d-15   ! [kg] : min mass, D_min=2um
  real(8), parameter :: xr_min = 2.60d-10   ! [kg] : min mass, D_min=79um
  real(8), parameter :: xi_min = 3.382d-13  ! [kg] : min mass, D_min=10um
  real(8), parameter :: xs_min = 1.847d-12  ! [kg] : min mass, D_min=20um
  real(8), parameter :: xg_min = 1.230d-10  ! [kg] : min mass, D_min=100um 
  ! refer to Seifert(2002) (Dr. Thesis, Table.5.1)
  real(8), parameter :: xc_max = 2.6d-10    ! [kg] : max, D_max=79um
  real(8), parameter :: xr_max = 5.00d-6    ! [kg] : max, D_max=2mm
  real(8), parameter :: xi_max = 1.377d-6   ! [kg] : max, D_max=5mm
  real(8), parameter :: xs_max = 7.519d-6   ! [kg] : max, D_max=1cm
  real(8), parameter :: xg_max = 4.900d-5   ! [kg] : max, D_max=1cm
  ! filter similar to Ikawa et al.(1991) sec.3.5
  real(8), parameter :: xmin_filter= xc_min
  ! filter of effective radius(1 micron)
  real(8), parameter :: rmin_re= 1.d-6
  !
  ! SB06(95),(96)
  real(8), parameter :: n0r_min= 2.5e5    ! [m-4]: min intercept parameter of rain
  real(8), parameter :: n0r_max= 2.0d7    ! [m-4]: max 
  real(8), parameter :: lambdar_min= 1.d3 ! [m-1]: min slope parameter of rain
  real(8), parameter :: lambdar_max= 1.d4 ! [m-1]: max 
  ! empirical value from Meyers etal.(1991), 1[/liter] = 1.d3[/m3]
  real(8), parameter :: nc_min = 1.d4     ! [m-3] empirical T.Mitsui
  real(8), parameter :: nr_min = 1.d0     ! [m-3] 1/1000 [/liter]
  real(8), parameter :: ni_min = 1.d0     ! [m-3]
  real(8), parameter :: ns_min = 1.d-4    ! [m-3] 
  real(8), parameter :: ng_min = 1.d-4    ! [m-3] 
  ! empirical filter 
  real(8), parameter :: lc_min = xc_min*nc_min
  real(8), parameter :: lr_min = xr_min*nr_min
  real(8), parameter :: li_min = xi_min*ni_min
  real(8), parameter :: ls_min = xs_min*ns_min
  real(8), parameter :: lg_min = xg_min*ng_min
  ! 
  real(8), parameter :: x_sep   = 2.6d-10 ! boundary mass between cloud and rain
  !
  real(8), parameter :: tem_min=100.D0 
  real(8), parameter :: rho_min=1.d-5     ! 3.e-3 is lower limit recognized in many experiments.
  real(8), parameter :: rhoi   = 916.7d0 
  ! for Seifert(2008)
  ! work parameter for gamma function, imported from MISC_gammafunc 
  real(8), parameter :: gfac_coef(6)=(/&
       +76.18009172947146D0, -86.50532032941677D0, +24.01409824083091D0,&
       -1.231739572450155D0, +0.1208650973866179D-2, -0.5395239384953D-5 /)
  real(8), parameter :: gfac_ser0=1.000000000190015D0
  !
  integer, save, private :: ntmax_phase_change = 1
  integer, save, private :: ntmax_collection   = 1
  integer, save, private :: ntmax_sedimentation= 1 ! 10/08/03 [Add] T.Mitsui
  !
  !--- standard density
  real(8), parameter :: rho_0 = 1.28D0
  !--- max number of Nc( activatable aerosol number concentration )
  real(8), allocatable, private, save :: nc_uplim(:,:)
  !
  !--- thermal conductivity of air
  real(8), parameter :: Ka0  = 2.428D-2 
  !<--- J/m/s/K : 0C/1atm
  real(8), parameter :: dKa_dT = 7.47D-5
  !<--- J/m/s/K/K : dependency of temperature
  !====== Ka = Ka0 + temc*dKa_dT
  !
  !--- Dynamic viscosity 
  real(8), parameter :: mua0 = 1.718D-5
  !<--- Kg/m/s : 0C/1atm
  real(8), parameter :: dmua_dT = 5.28D-8
  !<--- Kg/m/s/K : dependency of temperature
  !======  mua = mua0 + temc*dmua_dT
  !
  real(8), save, private :: xc_ccn = 1.d-12 ! [kg]
  real(8), save, private :: xi_ccn = 1.d-12 ! [kg] ! [move] 11/08/30 T.Mitsui
  !
  ! capacity of diffusional growth
  ! ( dependent of their geometries )
  real(8), save, private :: cap(HYDRO_MAX) 
  !
  ! constants for Diameter-Mass relation
  ! D = a * x^b
  real(8), save, private :: a_m(HYDRO_MAX) 
  real(8), save, private :: b_m(HYDRO_MAX) 
  ! constants for Terminal velocity-Mass relation
  ! vt = alpha * x^beta * f
  real(8), save, private :: alpha_v(HYDRO_MAX,2) 
  real(8), save, private :: beta_v(HYDRO_MAX,2) 
  real(8), save, private :: alpha_vn(HYDRO_MAX,2)   !
  real(8), save, private :: beta_vn(HYDRO_MAX,2)    !
  real(8), save, private :: gamma_v(HYDRO_MAX) 
  !  Aerodynamical factor for correction of terminal velocity.(Heymsfield and Iaquinta, 2000)
  !  vt(tem,pre) = vt0 * (pre/pre0)**a_pre0 * (tem/tem0)**a_tem0
  real(8), parameter :: pre0_vt   = 300.d2  ! 300hPa
  real(8), parameter :: tem0_vt   = 233.d0  ! -40degC
  real(8), parameter :: a_pre0_vt = -0.178d0
  real(8), parameter :: a_tem0_vt = -0.394d0
  ! Parameters to determine Droplet Size Distribution 
  ! as a General Gamma Distribution
  ! f(x) = A x^nu exp(-lambda x^mu )
  ! for Marshall Palmer Distribution ( popular for rain )
  !     mu=1/3, nu=-2/3
  ! for Gamma Distribution ( popular for cloud )
  !     mu=1
  real(8), save, private :: nu(HYDRO_MAX) 
  real(8), save, private :: mu(HYDRO_MAX) 
  ! Mitchell(1996), JAS, vol.53, No.12, pp.1710-1723
  !  area = a_area*D^b_area
  !  area = ax_area*x^bx_area
  ! Auer and Veal(1970), JAS, vol.27, pp.919-pp.926
  ! height = a_h*x^b_h( based on h=a_ar*D^b_ar,  ar:aspect ratio)
  real(8), save, private :: a_area(HYDRO_MAX)       !
  real(8), save, private :: b_area(HYDRO_MAX)       !
  real(8), save, private :: ax_area(HYDRO_MAX)      ! 
  real(8), save, private :: bx_area(HYDRO_MAX)      ! 
  ! parameters for radius of equivalent area
  ! r_ea = a_rea*x**b_rea
  real(8), save, private :: a_rea(HYDRO_MAX)        ! 
  real(8), save, private :: b_rea(HYDRO_MAX)        ! 
  real(8), save, private :: a_rea2(HYDRO_MAX)       ! 
  real(8), save, private :: b_rea2(HYDRO_MAX)       ! 
  real(8), save, private :: a_rea3(HYDRO_MAX)       ! 
  real(8), save, private :: b_rea3(HYDRO_MAX)       ! 
  !
  real(8), save, private :: a_d2vt(HYDRO_MAX)       ! 
  real(8), save, private :: b_d2vt(HYDRO_MAX)       ! 
  ! coefficient of x^2 moment of DSD
  ! Z = integral x*x*f(x) dx
  !   = coef_m2*N*(L/N)^2
  real(8), save, private :: coef_m2(HYDRO_MAX) 
  ! radar reflectivity coefficient defined by diameter
  real(8), save, private :: coef_d6(HYDRO_MAX)      !
  ! volume coefficient defined by diameter
  real(8), save, private :: coef_d3(HYDRO_MAX)      !
  ! coefficient of weighted mean diameter
  real(8), save, private :: coef_d(HYDRO_MAX) 
  ! coefficient of weighted mean d*d*v 
  real(8), save, private :: coef_d2v(HYDRO_MAX)     !
  ! coefficient of moment of d*d*v 
  real(8), save, private :: coef_md2v(HYDRO_MAX)    !
  !
  ! for effective radius(spherical particle)
  real(8), save, private :: coef_r2(HYDRO_MAX) 
  real(8), save, private :: coef_r3(HYDRO_MAX) 
  real(8), save, private :: coef_re(HYDRO_MAX) 
  ! for effective radius(hexagonal plate)
  real(8), save, private :: coef_rea2(HYDRO_MAX)    ! 
  real(8), save, private :: coef_rea3(HYDRO_MAX)    ! 
  logical, save, private :: opt_M96_ice=.true.               !
  logical, save, private :: opt_M96_column_ice=.false.       !
  !
  ! coefficeint of weighted mean terminal velocity
  ! vt0 is number weighted and
  ! vt1 is mass   weighted.
  real(8), save, private :: coef_vt0(HYDRO_MAX,2) 
  real(8), save, private :: coef_vt1(HYDRO_MAX,2) 
  real(8), save, private :: coef_deplc 
  real(8), save, private :: coef_dave_N(HYDRO_MAX) ! 
  real(8), save, private :: coef_dave_L(HYDRO_MAX) ! 
  ! diameter of terminal velocity branch
  !
  real(8), private, save :: d0_ni=261.76d-6   ! 
  real(8), private, save :: d0_li=398.54d-6   !
  real(8), parameter     :: d0_ns=270.03d-6   !
  real(8), parameter     :: d0_ls=397.47d-6   !
  real(8), parameter     :: d0_ng=269.08d-6   !
  real(8), parameter     :: d0_lg=376.36d-6   !
  !
  real(8), parameter :: coef_vtr_ar1=9.65d0    ! coef. for large branch
  ! original parameter of Rogers etal.(1993)
  real(8), parameter :: coef_vtr_br1=10.43d0    ! ...
  real(8), parameter :: coef_vtr_cr1=600.D0    ! ... 
  real(8), parameter :: coef_vtr_ar2=4.d3      ! coef. for small branch
  real(8), parameter :: coef_vtr_br2=12.d3     ! ...
  real(8), parameter :: d_vtr_branch=0.745d-3  ! 0.745 mm (diameter dividing 2-branches)
  ! equilibrium diameter of rain break-up
  real(8), parameter, private :: dr_eq   = 1.10d-3 ! eqilibrium diameter, Seifert 2008(36)
  ! coefficient of General Gamma.
  !  f(x)  = A x^nu exp(-lambda x^mu )
  ! lambda = coef_lambda * (L/N)^{-mu}  
  !     A  = coef_A*N*lambda^slope_A
  real(8), save, private :: coef_A(HYDRO_MAX) 
  real(8), save, private :: coef_lambda(HYDRO_MAX) 
  real(8), save, private :: slope_A(HYDRO_MAX)   
  ! coefficeint of weighted ventilation effect.
  ! large, and small branch is by PK97(13-60),(13-61),(13-88),(13-89) 
  real(8), save, private :: ah_vent  (HYDRO_MAX,2) ! 
  real(8), save, private :: bh_vent  (HYDRO_MAX,2) ! 
  real(8), save, private :: ah_vent0 (HYDRO_MAX,2) ! 
  real(8), save, private :: bh_vent0 (HYDRO_MAX,2) ! 
  real(8), save, private :: ah_vent1 (HYDRO_MAX,2) ! 
  real(8), save, private :: bh_vent1 (HYDRO_MAX,2) ! 
  ! coefficient of collision growth 
  real(8), save, private :: delta_b0 (HYDRO_MAX) 
  real(8), save, private :: delta_b1 (HYDRO_MAX) 
  real(8), save, private :: delta_ab0(HYDRO_MAX,HYDRO_MAX) 
  real(8), save, private :: delta_ab1(HYDRO_MAX,HYDRO_MAX) 
  !
  real(8), save, private :: theta_b0 (HYDRO_MAX) 
  real(8), save, private :: theta_b1 (HYDRO_MAX) 
  real(8), save, private :: theta_ab0(HYDRO_MAX,HYDRO_MAX) 
  real(8), save, private :: theta_ab1(HYDRO_MAX,HYDRO_MAX) 
  !
  logical, save, private :: opt_debug=.false.
  !
  logical, save, private :: opt_debug_tem=.false.
  logical, save, private :: opt_debug_inc=.true.
  logical, save, private :: opt_debug_act=.true.
  logical, save, private :: opt_debug_ree=.true.
  logical, save, private :: opt_debug_bcs=.true.

  integer, private, save :: MP_NSTEP_SEDIMENTATION
  real(8), private, save :: MP_RNSTEP_SEDIMENTATION
  real(8), private, save :: MP_DTSEC_SEDIMENTATION

  !
  ! metrics of vertical coordinate
  ! not used in SCALE-LES
  !
  real(8), private, allocatable, save :: gsgam2 (:,:)
  real(8), private, allocatable, save :: gsgam2h(:,:)
  real(8), private, allocatable, save :: gam2   (:,:)
  real(8), private, allocatable, save :: gam2h  (:,:)
  real(8), private, allocatable, save :: rgsgam2(:,:)
  real(8), private, allocatable, save :: rgs    (:,:)
  real(8), private, allocatable, save :: rgsh   (:,:)

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Setup Cloud Microphysics
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_setup
    use mod_time, only: &
       TIME_DTSEC_ATMOS_PHY_MP
    use mod_grid, only: &
       IMAX => GRID_IMAX, &
       JMAX => GRID_JMAX, &
       IJA  => GRID_IJA,  &
       KA   => GRID_KA
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Cloud Microphisics]/Categ[ATMOS]'
    if( IO_L ) write(IO_FID_LOG,*) '*** Wrapper for NDW6'

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

    call mp_ndw6_init( IMAX*JMAX )

    MP_NSTEP_SEDIMENTATION  = ntmax_sedimentation
    MP_RNSTEP_SEDIMENTATION = 1.D0 / real(ntmax_sedimentation,kind=8)
    MP_DTSEC_SEDIMENTATION  = TIME_DTSEC_ATMOS_PHY_MP * MP_RNSTEP_SEDIMENTATION

    allocate( gsgam2 (IJA,KA) )
    allocate( gsgam2h(IJA,KA) )
    allocate( gam2   (IJA,KA) )
    allocate( gam2h  (IJA,KA) )
    allocate( rgsgam2(IJA,KA) )
    allocate( rgs    (IJA,KA) )
    allocate( rgsh   (IJA,KA) )
    gsgam2 (:,:) = 1.D0
    gsgam2h(:,:) = 1.D0
    gam2   (:,:) = 1.D0
    gam2h  (:,:) = 1.D0
    rgsgam2(:,:) = 1.D0
    rgs    (:,:) = 1.D0
    rgsh   (:,:) = 1.D0

    return
  end subroutine ATMOS_PHY_MP_setup

  !-----------------------------------------------------------------------------
  !> Cloud Microphysics
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP
    use mod_grid, only : &
       KA   => GRID_KA,   &
       IA   => GRID_IA,   &
       JA   => GRID_JA,   &
       KS   => GRID_KS,   &
       KE   => GRID_KE
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait, &
       COMM_total
    use mod_atmos_vars, only: &
       var => atmos_var, &
       A_NAME,      &
       VA  => A_VA
    implicit none

    integer :: i, j, iv
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Microphysics'

    call mp_ndw6

    do iv = 1, VA
       ! fill IHALO & JHALO
       call COMM_vars8( var(:,:,:,iv), iv )
       call COMM_wait ( var(:,:,:,iv), iv )

       ! fill KHALO
       do j  = 1, JA
       do i  = 1, IA
          var(   1:KS-1,i,j,iv) = var(KS,i,j,iv)
          var(KE+1:KA,  i,j,iv) = var(KE,i,j,iv)
       enddo
       enddo
    enddo

    ! check total mass
    call COMM_total( var(:,:,:,:), A_NAME(:) )

    return
  end subroutine ATMOS_PHY_MP

  !-----------------------------------------------------------------------------
  subroutine mp_ndw6_init ( ijdim )
    use mod_stdio, only: &
       IO_FID_CONF,  &
       IO_FID_LOG,   &
       IO_L
    use mod_process, only: &
       PRC_MPIstop
    use mod_const, only : &
       CNST_PI    => CONST_PI, &
       CNST_UNDEF => CONST_UNDEF8
    use mod_atmos_vars, only: &
       NQW_STR => A_QWS, &
       NQW_END => A_QWE, &
       I_QV, &
       I_QC, &
       I_QR, &
       I_QI, &
       I_QS, &
       I_QG
    implicit none
    !
    integer, intent(in) :: ijdim
    !
    real(8), allocatable :: w1(:),w2(:),w3(:),w4(:),w5(:),w6(:),w7(:),w8(:)
    ! work for calculation of capacity, Mitchell and Arnott (1994) , eq.(9) 
    real(8) :: ar_ice_fix = 0.7d0   
    real(8) :: wcap1, wcap2
    ! work for ventilation coefficient
    logical :: flag_vent0(HYDRO_MAX), flag_vent1(HYDRO_MAX)
    integer :: ierr
    integer :: iw, ia, ib
    integer :: n
    !
    namelist /nm_mp_ndw6_init/       &
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
    namelist /nm_mp_ndw6_particles/ &
         a_m, b_m, alpha_v, beta_v, gamma_v, &
         alpha_vn, beta_vn,    &
         a_area, b_area, cap,  &
         nu, mu,               &
         opt_M96_column_ice,   &
         opt_M96_ice,          & 
         ar_ice_fix
    real(8), parameter :: eps_gamma=1.d-30

    a_m(:)         = CNST_UNDEF
    b_m(:)         = CNST_UNDEF
    alpha_v(:,:)   = CNST_UNDEF
    beta_v(:,:)    = CNST_UNDEF
    alpha_vn(:,:)  = CNST_UNDEF
    beta_vn(:,:)   = CNST_UNDEF
    gamma_v(:)     = CNST_UNDEF
    a_d2vt(:)      = CNST_UNDEF
    b_d2vt(:)      = CNST_UNDEF
    a_area(:)      = CNST_UNDEF
    b_area(:)      = CNST_UNDEF
    ax_area(:)     = CNST_UNDEF
    bx_area(:)     = CNST_UNDEF
    a_rea(:)       = CNST_UNDEF
    b_rea(:)       = CNST_UNDEF
    a_rea2(:)      = CNST_UNDEF
    b_rea2(:)      = CNST_UNDEF
    a_rea3(:)      = CNST_UNDEF
    b_rea3(:)      = CNST_UNDEF
    nu(:)          = CNST_UNDEF
    mu(:)          = CNST_UNDEF
    cap(:)         = CNST_UNDEF
    coef_m2(:)     = CNST_UNDEF
    coef_dave_N(:) = CNST_UNDEF
    coef_dave_L(:) = CNST_UNDEF
    coef_d(:)      = CNST_UNDEF
    coef_d3(:)     = CNST_UNDEF
    coef_d6(:)     = CNST_UNDEF
    coef_d2v(:)    = CNST_UNDEF
    coef_md2v(:)   = CNST_UNDEF
    coef_r2(:)     = CNST_UNDEF
    coef_r3(:)     = CNST_UNDEF
    coef_re(:)     = CNST_UNDEF
    coef_rea2(:)   = CNST_UNDEF
    coef_rea3(:)   = CNST_UNDEF
    coef_A(:)      = CNST_UNDEF
    slope_A(:)     = CNST_UNDEF
    coef_lambda(:) = CNST_UNDEF
    coef_vt0(:,:)  = CNST_UNDEF
    coef_vt1(:,:)  = CNST_UNDEF
    delta_b0(:)    = CNST_UNDEF
    delta_b1(:)    = CNST_UNDEF
    delta_ab0(:,:) = CNST_UNDEF
    delta_ab1(:,:) = CNST_UNDEF
    theta_b0(:)    = CNST_UNDEF
    theta_b1(:)    = CNST_UNDEF
    theta_ab0(:,:) = CNST_UNDEF
    theta_ab1(:,:) = CNST_UNDEF
    !
    ah_vent(:,:)   = CNST_UNDEF
    ah_vent0(:,:)  = CNST_UNDEF
    ah_vent1(:,:)  = CNST_UNDEF
    bh_vent(:,:)   = CNST_UNDEF
    bh_vent0(:,:)  = CNST_UNDEF
    bh_vent1(:,:)  = CNST_UNDEF
    !

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[NDW6]/Categ[ATMOS]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=nm_mp_ndw6_init,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist nm_mp_ndw6_init. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=nm_mp_ndw6_init)

    !
    ! default setting
    !
    ! Area parameters with mks unit originated by Mitchell(1996)
    a_area(I_QC) = CNST_PI/4.d0 ! sphere
    a_area(I_QR) = CNST_PI/4.d0 ! sphere
    a_area(I_QI) = 0.65d0*1.d-4*100.D0**(2.00d0)   ! Mitchell(1996), Hexagonal Plate
    a_area(I_QS) = 0.2285d0*1.d-4*100.D0**(1.88d0) ! Mitchell(1996), Aggregates
    a_area(I_QG) = 0.50d0*1.d-4*100.D0**(2.0d0)    ! Mitchell(1996), Lump Graupel
    b_area(I_QC) = 2.d0
    b_area(I_QR) = 2.d0
    b_area(I_QI) = 2.d0
    b_area(I_QS) = 1.88d0
    b_area(I_QG) = 2.d0
    !
    ! Seifert and Beheng(2006), Table. 1 or List of symbols
    !----------------------------------------------------------
    ! Diameter-Mass relationship
    ! D = a * x^b
    a_m(I_QC) = 0.124d0
    a_m(I_QR) = 0.124d0
    a_m(I_QI) = 0.217d0
    a_m(I_QS) = 8.156d0
    a_m(I_QG) = 0.190d0
    b_m(I_QC) = 1.d0/3.d0
    b_m(I_QR) = 1.d0/3.d0
    b_m(I_QI) = 0.302d0
    b_m(I_QS) = 0.526d0
    b_m(I_QG) = 0.323d0
    !----------------------------------------------------------
    ! Terminal velocity-Mass relationship
    ! vt = alpha * x^beta * (rho0/rho)^gamma
    alpha_v(I_QC,:)= 3.75d5
    alpha_v(I_QR,:)= 159.d0 ! not for sedimantation
    alpha_v(I_QI,:)= 317.d0
    alpha_v(I_QS,:)= 27.7d0
    alpha_v(I_QG,:)= 40.D0
    beta_v(I_QC,:) = 2.d0/3.d0
    beta_v(I_QR,:) = 0.266d0 ! not for sedimantation
    beta_v(I_QI,:) = 0.363d0
    beta_v(I_QS,:) = 0.216d0
    beta_v(I_QG,:) = 0.230d0
    gamma_v(I_QC)  = 1.d0
    ! This is high Reynolds number limit(Beard 1980)
    gamma_v(I_QR)  = 1.d0/2.d0 
    gamma_v(I_QI)  = 1.d0/2.d0 
    gamma_v(I_QS)  = 1.d0/2.d0 
    gamma_v(I_QG)  = 1.d0/2.d0 
    !----------------------------------------------------------       
    ! DSD parameters
    ! f(x) = A x^nu exp( -lambda x^mu )
    ! Gamma Disribution           : mu=1  , nu:arbitrary
    ! Marshall-Palmer Distribution: mu=1/3, nu:-2/3
    ! In the case of MP, f(D) dD = f(x)dx
    !                    f(x)    = c * f(D)/D^2 (c:coefficient)
    nu(I_QC) =  1.d0          ! arbitrary for Gamma
    nu(I_QR) = -1.d0/3.d0     ! nu(diameter)=1, equilibrium condition.
    nu(I_QI) =  1.d0          ! 
    nu(I_QS) =  1.d0          ! 
    nu(I_QG) =  1.d0          ! 
    !
    mu(I_QC) = 1.d0           ! Gamma
    mu(I_QR) = 1.d0/3.d0      ! Marshall Palmer
    mu(I_QI) = 1.d0/3.d0      ! 
    mu(I_QS) = 1.d0/3.d0      ! 
    mu(I_QG) = 1.d0/3.d0      ! 
    !----------------------------------------------------------
    ! Geomeries for diffusion growth
    ! Pruppacher and Klett(1997), (13-77)-(13-80) and
    ! originally derived by McDonald(1963b)
    ! sphere: cap=2
    ! plate : cap=pi
    ! needle with aspect ratio a/b
    !       : cap=log(2*a/b)
    cap(I_QC) = 2.d0    ! sphere
    cap(I_QR) = 2.d0    ! sphere
    cap(I_QI) = CNST_PI ! hexagonal plate
    cap(I_QS) = 2.d0    ! mix aggregates
    cap(I_QG) = 2.d0    ! lump
    !
    alpha_vn(:,:) = alpha_v(:,:)
    beta_vn(:,:) = beta_v(:,:)
    !------------------------------------------------------------------------
    !
    ! additional setting
    !

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=nm_mp_ndw6_particles,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist nm_mp_ndw6_particles. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=nm_mp_ndw6_particles)

    ! [Add] 10/08/03 T.Mitsui
    ! particles shapes are 
    if( opt_M96_ice ) then
       ! ice is randomly oriented Hexagonal plate (Auer and Veal 1970, Takano and Liou 1995, Mitchell 1996)
       ! snow is assemblages of planar polycrystals(Mitchell 1996)
       ! graupel is Lump graupel(R4b) is assumed(Mitchell 1996)
       a_area(I_QI)   = 0.120284936d0
       a_area(I_QS)   = 0.131488d0
       a_area(I_QG)   = 0.5d0
       b_area(I_QI)   = 1.850000d0
       b_area(I_QS)   = 1.880000d0
       b_area(I_QG)   = 2.d0     
       a_m(I_QI)      = 1.23655360084766d0
       a_m(I_QS)      = a_m(I_QI)
       a_m(I_QG)      = 0.346111225718402d0
       b_m(I_QI)      = 0.408329930583912d0 
       b_m(I_QS)      = b_m(I_QI)
       b_m(I_QG)      = 0.357142857142857d0
       !
       if( opt_M96_column_ice )then
          d0_ni=240.49d-6   ! this is column
          d0_li=330.09d-6   ! this is column
          a_area(I_QI)= (0.684d0*1.d-4)*10.D0**(2.d0*2.0d0)
          b_area(I_QI)= 2.0d0
          a_m(I_QI)   = 0.19834046116844d0
          b_m(I_QI)   = 0.343642611683849d0
          ! [Add] 11/08/30 T.Mitsui
          ! approximated by the capacity for prolate spheroid with constant aspect ratio
          wcap1       = sqrt(1.d0-ar_ice_fix**2)
          wcap2       = log( (1.d0+wcap1)/ar_ice_fix )
          cap(I_QI)   = 2.d0*wcap2/wcap1
          !
       end if
       !
       ! These value are derived by least-square fitting in the range
       ! qi [100um:1000um] in diameter
       ! qs [100um:1000um] in diameter
       ! qg [200um:2000um] in diameter
       !                      small branch          , large branch
       alpha_v (I_QI,:) =(/ 5798.60107421875d0, 167.347076416016d0/)
       alpha_vn(I_QI,:) =(/ 12408.1777343750d0, 421.799865722656d0/)          
       if( opt_M96_column_ice )then
          alpha_v (I_QI,:) = (/2901.0d0, 32.20d0/)
          alpha_vn(I_QI,:) = (/9675.2d0, 64.16d0/)
       end if
       alpha_v (I_QS,:)    =(/ 15173.3916015625d0, 305.678619384766d0/)
       alpha_vn(I_QS,:)    =(/ 29257.1601562500d0, 817.985717773438d0/)
       alpha_v (I_QG,:)    =(/ 15481.6904296875d0, 311.642242431641d0/)
       alpha_vn(I_QG,:)    =(/ 27574.6562500000d0, 697.536132812500d0/)
       !
       beta_v (I_QI,:)  =(/ 0.504873454570770d0, 0.324817866086960d0/)
       beta_vn(I_QI,:)  =(/ 0.548495233058929d0, 0.385287821292877d0/)
       if( opt_M96_column_ice )then
          beta_v (I_QI,:)  =(/ 0.465552181005478d0, 0.223826110363007d0/)
          beta_vn(I_QI,:)  =(/ 0.530453503131866d0, 0.273761242628098d0/)          
       end if
       beta_v (I_QS,:)     =(/ 0.528109610080719d0, 0.329863965511322d0/)
       beta_vn(I_QS,:)     =(/ 0.567154467105865d0, 0.393876969814301d0/)
       beta_v (I_QG,:)     =(/ 0.534656763076782d0, 0.330253750085831d0/)
       beta_vn(I_QG,:)     =(/ 0.570551633834839d0, 0.387124240398407d0/)       
    end if
    ! 
    ! area-diameter relation => area-mass relation
    ax_area(I_QC:I_QG) = a_area(I_QC:I_QG)*a_m(I_QC:I_QG)**b_area(I_QC:I_QG)
    bx_area(I_QC:I_QG) = b_area(I_QC:I_QG)*b_m(I_QC:I_QG)
    ! 
    ! radius of equivalent area - m ass relation
    ! pi*rea**2 = ax_area*x**bx_area
    a_rea(I_QC:I_QG)   = sqrt(ax_area(I_QC:I_QG)/CNST_PI)
    b_rea(I_QC:I_QG)   = bx_area(I_QC:I_QG)/2.d0
    a_rea2(I_QC:I_QG)  = a_rea(I_QC:I_QG)**2
    b_rea2(I_QC:I_QG)  = b_rea(I_QC:I_QG)*2.d0
    a_rea3(I_QC:I_QG)  = a_rea(I_QC:I_QG)**3
    b_rea3(I_QC:I_QG)  = b_rea(I_QC:I_QG)*3.d0
    !
    a_d2vt(I_QC:I_QG)=alpha_v(I_QC:I_QG,2)*(1.d0/alpha_v(I_QC:I_QG,2))**(beta_v(I_QC:I_QG,2)/b_m(I_QC:I_QG))
    b_d2vt(I_QC:I_QG)=(beta_v(I_QC:I_QG,2)/b_m(I_QC:I_QG))
    !
    ! Calculation of Moment Coefficient
    !
    allocate( w1(I_QC:I_QG), w2(I_QC:I_QG), w3(I_QC:I_QG), w4(I_QC:I_QG) )
    allocate( w5(I_QC:I_QG), w6(I_QC:I_QG), w7(I_QC:I_QG), w8(I_QC:I_QG) )
    w1(:) = 0.D0
    w2(:) = 0.D0
    w3(:) = 0.D0
    w4(:) = 0.D0
    w5(:) = 0.D0
    w6(:) = 0.D0
    w7(:) = 0.D0
    w8(:) = 0.D0
    !-------------------------------------------------------
    ! moment coefficient
    ! SB06 (82)
    ! M^n = coef_mn * N * (L/N)**n       
    ! M^2 = Z = coef_m2 * N *(L/N)**2
    ! a*M^b = a*integral x^b f(x) dx = ave D
    do iw=I_QC,I_QG
       n = 2
       w1(iw) = gammafunc( (n+nu(iw)+1.d0)/mu(iw) )
       w2(iw) = gammafunc( (nu(iw)+1.d0)/mu(iw) ) 
       w3(iw) = gammafunc( (nu(iw)+2.d0)/mu(iw) )
       coef_m2(iw) = w1(iw)/w2(iw)*(  w2(iw)/w3(iw)  )**n
       ! 
       w4(iw) = gammafunc( (b_m(iw)+nu(iw)+1.d0)/mu(iw) )
       coef_d(iw) = a_m(iw) * w4(iw)/w2(iw)*(  w2(iw)/w3(iw)  )**b_m(iw)
       w5(iw) = gammafunc( (2.d0*b_m(iw)+beta_v(iw,2)+nu(iw)+1.d0)/mu(iw) )
       w6(iw) = gammafunc( (3.d0*b_m(iw)+beta_v(iw,2)+nu(iw)+1.d0)/mu(iw) )
       coef_d2v(iw) = a_m(iw) * w6(iw)/w5(iw)* ( w2(iw)/w3(iw) )**b_m(iw)
       coef_md2v(iw)=           w5(iw)/w2(iw)* ( w2(iw)/w3(iw) )**(2.d0*b_m(iw)+beta_v(iw,2))       
       ! 09/04/14 [Add] T.Mitsui, volume and radar reflectivity
       w7(iw) = gammafunc( (3.d0*b_m(iw)+nu(iw)+1.d0)/mu(iw) )
       coef_d3(iw)  = a_m(iw)**3 * w7(iw)/w2(iw)*(  w2(iw)/w3(iw)  )**(3.d0*b_m(iw))
       w8(iw) = gammafunc( (6.d0*b_m(iw)+nu(iw)+1.d0)/mu(iw) )
       coef_d6(iw)  = a_m(iw)**6 * w8(iw)/w2(iw)*(  w2(iw)/w3(iw)  )**(6.d0*b_m(iw))
    end do
    !
    coef_deplc = coef_d(I_QC)/a_m(I_QC)
    !-------------------------------------------------------
    ! coefficient of 2nd and 3rd moments for effective radius
    ! for spherical particle
    do iw=I_QC,I_QG
       ! integ r^2 f(x)dx
       w1(iw) = gammafunc( (2.d0*b_m(iw)+nu(iw)+1.d0)/mu(iw) )
       w2(iw) = gammafunc( (nu(iw)+1.d0)/mu(iw) ) 
       w3(iw) = gammafunc( (nu(iw)+2.d0)/mu(iw) )
       ! integ r^3 f(x)dx
       w4(iw) = gammafunc( (3.d0*b_m(iw)+nu(iw)+1.d0)/mu(iw) )
       !
       coef_r2(iw)=w1(iw)/w2(iw)*( w2(iw)/w3(iw) )**(2.d0*b_m(iw))
       coef_r3(iw)=w4(iw)/w2(iw)*( w2(iw)/w3(iw) )**(3.d0*b_m(iw))
       coef_re(iw)=coef_r3(iw)/coef_r2(iw)
       !
    end do
    !-------------------------------------------------------
    ! coefficient for effective radius of equivalent area and   
    ! coefficient for volume of equivalent area
    do iw=I_QC,I_QG
       w1(iw) = gammafunc( (nu(iw)+1.d0)/mu(iw) ) 
       w2(iw) = gammafunc( (nu(iw)+2.d0)/mu(iw) )
       w3(iw) = gammafunc( (b_rea2(iw)+nu(iw)+1.d0)/mu(iw) )
       w4(iw) = gammafunc( (b_rea3(iw)+nu(iw)+1.d0)/mu(iw) ) 
       !
       coef_rea2(iw) = w3(iw)/w1(iw)*( w1(iw)/w2(iw) )**b_rea2(iw)
       coef_rea3(iw) = w4(iw)/w1(iw)*( w1(iw)/w2(iw) )**b_rea3(iw)
    end do
    !-------------------------------------------------------
    ! coefficient of gamma-distribution
    ! SB06(80)
    do iw=I_QC,I_QG
       w1(iw) = gammafunc( (nu(iw) + 1.d0)/mu(iw) )    
       w2(iw) = gammafunc( (nu(iw) + 2.d0)/mu(iw) )
       coef_A(iw)      = mu(iw)/w1(iw)
       slope_A(iw)     = w1(iw)
       coef_lambda(iw) = (w1(iw)/w2(iw))**(-mu(iw))
    end do
    !-------------------------------------------------------
    ! coefficient for terminal velocity in sedimentation
    ! SB06(78)
    do ia=1,2
       do iw=I_QC,I_QG
          n = 0
          w1(iw) = gammafunc( (beta_vn(iw,ia) + nu(iw) + 1.d0 + n)/mu(iw) )
          w2(iw) = gammafunc( (                 nu(iw) + 1.d0 + n)/mu(iw) )
          w3(iw) = gammafunc( (nu(iw) + 1.d0)/mu(iw) )
          w4(iw) = gammafunc( (nu(iw) + 2.d0)/mu(iw) )
          ! coefficient of terminal velocity for number
          coef_vt0(iw,ia)=alpha_vn(iw,ia)*w1(iw)/w2(iw)*(w3(iw)/w4(iw))**beta_vn(iw,ia)
          n = 1
          w1(iw) = gammafunc( (beta_v(iw,ia) + nu(iw) + 1.d0 + n)/mu(iw) )
          w2(iw) = gammafunc( (                nu(iw) + 1.d0 + n)/mu(iw) )
          ! coefficient of terminal velocity for mass
          coef_vt1(iw,ia)=alpha_v(iw,ia)*w1(iw)/w2(iw)*(w3(iw)/w4(iw))**beta_v(iw,ia)
       end do
    end do
    ! coefficient for weighted diameter used in calculation of terminal velocity
    do iw=I_QC,I_QG
       w1(iw) = gammafunc( (       b_m(iw) + nu(iw) + 1.d0)/mu(iw) )
       w2(iw) = gammafunc( (1.d0 + b_m(iw) + nu(iw) + 1.d0)/mu(iw) )
       w3(iw) = gammafunc( (nu(iw) + 1.d0)/mu(iw) )
       w4(iw) = gammafunc( (nu(iw) + 2.d0)/mu(iw) )       
       coef_dave_N(iw) =  (w1(iw)/w3(iw))*(w3(iw)/w4(iw))**      b_m(iw)         
       coef_dave_L(iw) =  (w2(iw)/w3(iw))*(w3(iw)/w4(iw))**(1.d0+b_m(iw)) 
    end do
    !-------------------------------------------------------
    !
    ah_vent(I_QC,1:2) = (/1.000d0,1.000d0/) ! no effect
    ah_vent(I_QR,1:2) = (/1.000d0,0.780d0/)
    ah_vent(I_QI,1:2) = (/1.000d0,0.860d0/)
    ah_vent(I_QS,1:2) = (/1.000d0,0.780d0/)
    ah_vent(I_QG,1:2) = (/1.000d0,0.780d0/)
    bh_vent(I_QC,1:2) = (/0.000d0,0.000d0/)
    bh_vent(I_QR,1:2) = (/0.108d0,0.308d0/)
    bh_vent(I_QI,1:2) = (/0.140d0,0.280d0/)
    bh_vent(I_QS,1:2) = (/0.108d0,0.308d0/)
    bh_vent(I_QG,1:2) = (/0.108d0,0.308d0/)    
    !
    do iw=I_QC,I_QG
       n = 0
       if( (nu(iw) + b_m(iw) + n) > eps_gamma  )then
          w1(iw) = gammafunc( (nu(iw) + b_m(iw) + n)/mu(iw) )
          w2(iw) = gammafunc( (nu(iw) + 1.d0)/mu(iw) )
          w3(iw) = gammafunc( (nu(iw) + 2.d0)/mu(iw) )    
          ah_vent0(iw,1)= ah_vent(iw,1)*(w1(iw)/w2(iw))*(w2(iw)/w3(iw))**(b_m(iw)+n-1.d0)
          ah_vent0(iw,2)= ah_vent(iw,2)*(w1(iw)/w2(iw))*(w2(iw)/w3(iw))**(b_m(iw)+n-1.d0)
          flag_vent0(iw)=.true.
       else
          ah_vent0(iw,1)= 1.d0
          ah_vent0(iw,2)= 1.d0
          flag_vent0(iw)=.false.
       end if
       n = 1
       if( (nu(iw) + b_m(iw) + n) > eps_gamma  )then
          w1(iw) = gammafunc( (nu(iw) + b_m(iw) + n)/mu(iw) )    
          w2(iw) = gammafunc( (nu(iw) + 1.d0)/mu(iw) )
          w3(iw) = gammafunc( (nu(iw) + 2.d0)/mu(iw) )    
          ah_vent1(iw,1)= ah_vent(iw,1)*(w1(iw)/w2(iw))*(w2(iw)/w3(iw))**(b_m(iw)+n-1.d0)
          ah_vent1(iw,2)= ah_vent(iw,2)*(w1(iw)/w2(iw))*(w2(iw)/w3(iw))**(b_m(iw)+n-1.d0)
          flag_vent1(iw)=.true.
       else
          ah_vent1(iw,1)= 1.d0
          ah_vent1(iw,2)= 1.d0
          flag_vent1(iw)=.true.
       end if
    end do
    do iw=I_QC,I_QG
       n = 0
       if( (nu(iw) + 1.5d0*b_m(iw) + 0.5d0*beta_v(iw,1) + n) < eps_gamma )then
          flag_vent0(iw)=.false.
       end if
       if(flag_vent0(iw))then
          w1(iw) = gammafunc( (nu(iw) + 1.5d0*b_m(iw) + 0.5d0*beta_v(iw,1) + n)/mu(iw) )
          w2(iw) = gammafunc( (nu(iw) + 1.d0)/mu(iw) )
          w3(iw) = gammafunc( (nu(iw) + 2.d0)/mu(iw) )    
          ! [Add] 11/08/30 T.Mitsui
          w4(iw) = gammafunc( (nu(iw) + 2.d0*b_m(iw) + beta_v(iw,1) + n)/mu(iw) )
          bh_vent0(iw,1)=bh_vent(iw,1)*(w4(iw)/w2(iw))*(w2(iw)/w3(iw))**(2.0d0*b_m(iw)+beta_v(iw,1)+n-1.d0)
          w5(iw) = gammafunc( (nu(iw) + 1.5d0*b_m(iw) + 0.5d0*beta_v(iw,2) + n)/mu(iw) )
          bh_vent0(iw,2)=bh_vent(iw,2)*(w5(iw)/w2(iw))*(w2(iw)/w3(iw))**(1.5d0*b_m(iw)+0.5d0*beta_v(iw,2)+n-1.d0)
       else
          bh_vent0(iw,1) = 0.D0
          bh_vent0(iw,2) = 0.D0
       end if
       !
       n = 1
       if( (nu(iw) + 1.5d0*b_m(iw) + 0.5d0*beta_v(iw,1) + n) < eps_gamma )then
          flag_vent1(iw)=.false.
       end if
       if(flag_vent1(iw))then
          w1(iw) = gammafunc( (nu(iw) + 1.5d0*b_m(iw) + 0.5d0*beta_v(iw,1) + n)/mu(iw) )
          w2(iw) = gammafunc( (nu(iw) + 1.d0)/mu(iw) )
          w3(iw) = gammafunc( (nu(iw) + 2.d0)/mu(iw) )    
          ! [Add] 11/08/30 T.Mitsui
          w4(iw) = gammafunc( (nu(iw) + 2.d0*b_m(iw) + beta_v(iw,1) + n)/mu(iw) )
          bh_vent1(iw,1)=bh_vent(iw,1)*(w4(iw)/w2(iw))*(w2(iw)/w3(iw))**(2.0d0*b_m(iw)+beta_v(iw,1)+n-1.d0)
          !
          w5(iw) = gammafunc( (nu(iw) + 1.5d0*b_m(iw) + 0.5d0*beta_v(iw,2) + n)/mu(iw) )
          bh_vent1(iw,2)=bh_vent(iw,2)*(w5(iw)/w2(iw))*(w2(iw)/w3(iw))**(1.5d0*b_m(iw)+0.5d0*beta_v(iw,2)+n-1.d0)
       else
          bh_vent1(iw,1) = 0.D0
          bh_vent1(iw,2) = 0.D0
       end if
    end do
    !-------------------------------------------------------
    ! coefficient for collision process
    ! stochastic coefficient for collision cross section
    ! sb06 (90) -- self collection
    do iw=I_QC,I_QG
       n = 0
       w1(iw) = gammafunc( (2.d0*b_rea(iw) + nu(iw) + 1.d0 + n)/mu(iw) )
       w2(iw) = gammafunc( (nu(iw) + 1.d0)/mu(iw) )
       w3(iw) = gammafunc( (nu(iw) + 2.d0)/mu(iw) )
       delta_b0(iw) = w1(iw)/w2(iw) &
            *( w2(iw)/w3(iw) )**(2.d0*b_rea(iw) + n) 
       n = 1
       w1(iw) = gammafunc( (2.d0*b_rea(iw) + nu(iw) + 1.d0 + n)/mu(iw) )
       delta_b1(iw) = w1(iw)/w2(iw) &
            *( w2(iw)/w3(iw) )**(2.d0*b_rea(iw) + n) 
    end do
    ! stochastic coefficient for collision cross section
    ! sb06(91) -- riming( collide with others )
    do iw=I_QC,I_QG
       n = 0
       w1(iw) = gammafunc( (b_rea(iw) + nu(iw) + 1.d0 + n)/mu(iw) )
       w2(iw) = gammafunc( (nu(iw) + 1.d0)/mu(iw) )
       w3(iw) = gammafunc( (nu(iw) + 2.d0)/mu(iw) )
       w4(iw) = gammafunc( (b_rea(iw) + nu(iw) + 1.d0    )/mu(iw) )
       n = 1
       w5(iw) = gammafunc( (b_rea(iw) + nu(iw) + 1.d0 + n)/mu(iw) )
    end do
    ! ia > ib ( larger particles "a" catch smaller particles "b" )
    do ia=I_QC,I_QG
       do ib=I_QC,I_QG
          n=0 ! 
          ! NOTE, collected  particle has a moment of n.
          !       collecting particle has only number(n=0).
          delta_ab0(ia,ib) = 2.d0*(w1(ib)/w2(ib))*(w4(ia)/w2(ia)) &
               * ( w2(ib)/w3(ib) )**(b_rea(ib)+n) &
               * ( w2(ia)/w3(ia) )**(b_rea(ia)  )
          n=1 ! 
          delta_ab1(ia,ib) = 2.d0*(w5(ib)/w2(ib))*(w4(ia)/w2(ia)) &
               * ( w2(ib)/w3(ib) )**(b_rea(ib)+n) &
               * ( w2(ia)/w3(ia) )**(b_rea(ia)  )
       end do
    end do
    ! stochastic coefficient for terminal velocity
    ! sb06(92) -- self collection
    ! assuming equivalent area circle.
    do iw=I_QC,I_QG
       n = 0
       w1(iw) = gammafunc( (2.d0*beta_v(iw,2) + 2.d0*b_rea(iw) + nu(iw) + 1.d0 + n)/mu(iw) )
       w2(iw) = gammafunc( (                        2.d0*b_rea(iw) + nu(iw) + 1.d0 + n)/mu(iw) )
       w3(iw) = gammafunc( (nu(iw) + 1.d0)/mu(iw) )
       w4(iw) = gammafunc( (nu(iw) + 2.d0)/mu(iw) )
       theta_b0(iw) = w1(iw)/w2(iw) * ( w3(iw)/w4(iw) )**(2.d0*beta_v(iw,2))
       n = 1
       w1(iw) = gammafunc( (2.d0*beta_v(iw,2) + 2.d0*b_rea(iw) + nu(iw) + 1.d0 + n)/mu(iw) )
       w2(iw) = gammafunc( (                        2.d0*b_rea(iw) + nu(iw) + 1.d0 + n)/mu(iw) )
       theta_b1(iw) = w1(iw)/w2(iw) * ( w3(iw)/w4(iw) )**(2.d0*beta_v(iw,2))
    end do
    !
    ! stochastic coefficient for terminal velocity
    ! sb06(93) -- riming( collide with others )
    do iw=I_QC,I_QG
       n = 0
       w1(iw) = gammafunc( (beta_v(iw,2) + 2.d0*b_rea(iw) + nu(iw) + 1.d0 + n)/mu(iw) )
       w2(iw) = gammafunc( (                   2.d0*b_rea(iw) + nu(iw) + 1.d0 + n)/mu(iw) )
       w3(iw) = gammafunc( (beta_v(iw,2) + 2.d0*b_rea(iw) + nu(iw) + 1.d0    )/mu(iw) )
       w4(iw) = gammafunc( (                   2.d0*b_rea(iw) + nu(iw) + 1.d0    )/mu(iw) )
       !
       w5(iw) = gammafunc( (nu(iw) + 1.d0)/mu(iw) )
       w6(iw) = gammafunc( (nu(iw) + 2.d0)/mu(iw) )
       n = 1
       w7(iw) = gammafunc( (beta_v(iw,2) + b_rea(iw) + nu(iw) + 1.d0 + n)/mu(iw) )
       w8(iw) = gammafunc( (                   b_rea(iw) + nu(iw) + 1.d0 + n)/mu(iw) )
    end do
    ! ia > ib ( larger particles "a" catch smaller particles "b" )
    do ia=I_QC,I_QG
       do ib=I_QC,I_QG
          theta_ab0(ia,ib) = 2.d0 * (w1(ib)/w2(ib))*(w3(ia)/w4(ia)) &
               * (w5(ia)/w6(ia))**beta_v(ia,2) &
               * (w5(ib)/w6(ib))**beta_v(ib,2) 
          theta_ab1(ia,ib) = 2.d0 * (w7(ib)/w8(ib))*(w3(ia)/w4(ia)) &
               * (w5(ia)/w6(ia))**beta_v(ia,2) &
               * (w5(ib)/w6(ib))**beta_v(ib,2) 
       end do
    end do
    !
    deallocate(w1,w2,w3,w4,w5,w6,w7,w8)
    !
    write(IO_FID_LOG,'(100a16)')     "LABEL       ",WLABEL(:)
    write(IO_FID_LOG,'(a,100e16.6)') "capacity    ",cap(:) ! [Add] 11/08/30 T.Mitsui
    write(IO_FID_LOG,'(a,100e16.6)') "coef_m2     ",coef_m2(:)
    write(IO_FID_LOG,'(a,100e16.6)') "coef_d      ",coef_d(:)
    !
    write(IO_FID_LOG,'(a,100e16.6)') "coef_d3     ",coef_d3(:)
    write(IO_FID_LOG,'(a,100e16.6)') "coef_d6     ",coef_d6(:)
    write(IO_FID_LOG,'(a,100e16.6)') "coef_d2v    ",coef_d2v(:)
    write(IO_FID_LOG,'(a,100e16.6)') "coef_md2v   ",coef_md2v(:)    
    write(IO_FID_LOG,'(a,100e16.6)') "a_d2vt      ",a_d2vt(:) 
    write(IO_FID_LOG,'(a,100e16.6)') "b_d2vt      ",b_d2vt(:) 
    !    
    write(IO_FID_LOG,'(a,100e16.6)') "coef_r2     ",coef_r2(:)
    write(IO_FID_LOG,'(a,100e16.6)') "coef_r3     ",coef_r3(:)
    write(IO_FID_LOG,'(a,100e16.6)') "coef_re     ",coef_re(:)
    !
    write(IO_FID_LOG,'(a,100e16.6)') "a_area      ",a_area(:)
    write(IO_FID_LOG,'(a,100e16.6)') "b_area      ",b_area(:)
    write(IO_FID_LOG,'(a,100e16.6)') "ax_area     ",ax_area(:)
    write(IO_FID_LOG,'(a,100e16.6)') "bx_area     ",bx_area(:)
    write(IO_FID_LOG,'(a,100e16.6)') "a_rea       ",a_rea(:)
    write(IO_FID_LOG,'(a,100e16.6)') "b_rea       ",b_rea(:)
    write(IO_FID_LOG,'(a,100e16.6)') "a_rea3      ",a_rea3(:)
    write(IO_FID_LOG,'(a,100e16.6)') "b_rea3      ",b_rea3(:)
    !
    write(IO_FID_LOG,'(a,100e16.6)') "coef_rea2   ",coef_rea2(:)
    write(IO_FID_LOG,'(a,100e16.6)') "coef_rea3   ",coef_rea3(:)
    write(IO_FID_LOG,'(a,100e16.6)') "coef_vt0    ",coef_vt0(:,1)
    write(IO_FID_LOG,'(a,100e16.6)') "coef_vt1    ",coef_vt1(:,1)
    write(IO_FID_LOG,'(a,100e16.6)') "coef_A      ",coef_A(:)
    write(IO_FID_LOG,'(a,100e16.6)') "coef_lambda ",coef_lambda(:)
    
    write(IO_FID_LOG,'(a,100e16.6)') "ah_vent0 sml",ah_vent0(:,1)
    write(IO_FID_LOG,'(a,100e16.6)') "ah_vent0 lrg",ah_vent0(:,2)
    write(IO_FID_LOG,'(a,100e16.6)') "ah_vent1 sml",ah_vent1(:,1)
    write(IO_FID_LOG,'(a,100e16.6)') "ah_vent1 lrg",ah_vent1(:,2)
    write(IO_FID_LOG,'(a,100e16.6)') "bh_vent0 sml",bh_vent0(:,1)
    write(IO_FID_LOG,'(a,100e16.6)') "bh_vent0 lrg",bh_vent0(:,2)
    write(IO_FID_LOG,'(a,100e16.6)') "bh_vent1 sml",bh_vent1(:,1)
    write(IO_FID_LOG,'(a,100e16.6)') "bh_vent1 lrg",bh_vent1(:,2)
    
    write(IO_FID_LOG,'(a,100e16.6)') "delta_b0    ",delta_b0(:)
    write(IO_FID_LOG,'(a,100e16.6)') "delta_b1    ",delta_b1(:)
    write(IO_FID_LOG,'(a,100e16.6)') "theta_b0    ",theta_b0(:)
    write(IO_FID_LOG,'(a,100e16.6)') "theta_b1    ",theta_b1(:)

    do ia=NQW_STR,NQW_END
       write(IO_FID_LOG,'(a,a10,a,100e16.6)') "delta0(a,b)=(",trim(WLABEL(ia)),",b)=",delta_ab0(ia,:)
    enddo
    do ia=NQW_STR,NQW_END
       write(IO_FID_LOG,'(a,a10,a,100e16.6)') "delta1(a,b)=(",trim(WLABEL(ia)),",b)=",delta_ab1(ia,:)
    enddo
    do ia=NQW_STR,NQW_END
       write(IO_FID_LOG,'(a,a10,a,100e16.6)') "theta0(a,b)=(",trim(WLABEL(ia)),",b)=",theta_ab0(ia,:)
    enddo
    do ia=NQW_STR,NQW_END
       write(IO_FID_LOG,'(a,a10,a,100e16.6)') "theta1(a,b)=(",trim(WLABEL(ia)),",b)=",theta_ab1(ia,:)
    enddo

    allocate(nc_uplim(ijdim,1))
    nc_uplim(:,:) = 150.d6

    return
  end subroutine mp_ndw6_init
  !-----------------------------------------------------------------------------
  subroutine mp_ndw6
    use mod_process, only: &
       l_region => PRC_myrank
    use mod_const, only : &
       RovCP => CONST_RovCP, &
       PRE00 => CONST_PRE00
    use mod_time, only: &
       dt => TIME_DTSEC_ATMOS_PHY_MP, &
       ct => TIME_NOWSEC
    use mod_grid, only: &
       kdim   => GRID_KA,   &
       IA   => GRID_IA,   &
       JA   => GRID_JA,   &
       IMAX => GRID_IMAX,   &
       JMAX => GRID_JMAX,   &
       kmin   => GRID_KS,   &
       kmax   => GRID_KE,   &
       IS   => GRID_IS,   &
       IE   => GRID_IE,   &
       JS   => GRID_JS,   &
       JE   => GRID_JE,   &
       z    => GRID_CZ,   &
       dz   => GRID_CDZ,   &
       ijdim  => GRID_IJA,  &
       IJS  => GRID_IJS,  &
       IJE  => GRID_IJE
    use mod_precipitation, only : &
       precipitation
    use mod_thrmdyn, only : &
       thrmdyn_cv, &
       thrmdyn_qd, &
       thrmdyn_cp, &
       thrmdyn_tempre, &
       thrmdyn_tempre2
    use mod_satadjust, only : &
       moist_psat_ice,     &
       moist_psat_water,   &
       moist_dqsw_dtem_rho
    use mod_atmos_vars, only: &
       var => atmos_var, &
       nqmax => A_QA,  &
       CVw => A_CVw, &
       I_DENS,      &
       I_MOMX,      &
       I_MOMY,      &
       I_MOMZ,      &
       I_RHOT,      &
       I_QV,        &
       I_QC,        &
       I_QR,        &
       I_QI,        &
       I_QS,        &
       I_QG,        &
       I_NC,        &
       I_NR,        &
       I_NI,        &
       I_NS,        &
       I_NG
    use mod_atmos_cnst, only :&
         CNST_CV,       &
         CNST_CVV,      &
         CNST_CL,       &
         CNST_CI,       &
         CNST_LH0,      &
         CNST_LHS0,     &
         CNST_LHF0,     &
         CNST_LH00,     &
         CNST_LHS00,    &
         CNST_LHF00,    &
         CNST_EGRAV,    &
         CNST_RVAP,     &
         CNST_RAIR,     & 
         CNST_EPSV,     &
         CNST_PSAT0,    &
         CNST_TEM00,    &
         CNST_PRES0,    &
         CNST_TMELT,    &  
         CNST_UNDEF,    &
         LHV,LHF,          & 
         NNW_STR, NNW_END, &
         NQW_STR, NQW_END
    implicit none

    !
    ! primary variables
    !
    real(8) :: rhogvx(ijdim,kdim)
    real(8) :: rhogvy(ijdim,kdim)
    real(8) :: rhogw (ijdim,kdim)
    real(8) :: rhogq (ijdim,kdim,nqmax)
    real(8) :: rhog  (ijdim,kdim)
    real(8) :: th    (ijdim,kdim)

    !
    ! diagnostic variables
    !
    real(8) :: qd(ijdim,kdim)
    real(8) :: q(ijdim,kdim,nqmax)
    real(8) :: rho(ijdim,kdim)
    real(8) :: tem(ijdim,kdim)
    real(8) :: pre(ijdim,kdim)
    real(8) :: rhoge(ijdim,kdim)
    real(8) :: w(ijdim,kdim)
    !
    real(8) :: lv(ijdim,kdim)
    real(8) :: lc(ijdim,kdim), nc(ijdim,kdim)
    real(8) :: lr(ijdim,kdim), nr(ijdim,kdim)
    real(8) :: li(ijdim,kdim), ni(ijdim,kdim)
    real(8) :: ls(ijdim,kdim), ns(ijdim,kdim)
    real(8) :: lg(ijdim,kdim), ng(ijdim,kdim)
    !
    real(8) :: xc(ijdim,kdim) ! lc/nc
    real(8) :: xr(ijdim,kdim) ! lr/nr
    real(8) :: xi(ijdim,kdim) ! li/ni
    real(8) :: xs(ijdim,kdim) ! ls/ns
    real(8) :: xg(ijdim,kdim) ! lg/ng
    !
    real(8) :: dc_xa(ijdim,kdim) ! 
    real(8) :: dr_xa(ijdim,kdim) ! 
    real(8) :: di_xa(ijdim,kdim) ! 
    real(8) :: ds_xa(ijdim,kdim) ! 
    real(8) :: dg_xa(ijdim,kdim) ! 
    real(8) :: vt_xa(ijdim,kdim,HYDRO_MAX,2) ! 


    real(8) :: wtem(ijdim,kdim)      ! filtered temperature 
    real(8) :: esw(ijdim,kdim)       ! saturated vapor pressure(water)
    real(8) :: esi(ijdim,kdim)       ! saturated vapor pressure(ice)
    real(8) :: r_cva                  ! work cva    
    !
    real(8) :: rho_fac               ! 
    real(8) :: rho_fac_c(ijdim,kdim) ! factor for cloud
    real(8) :: rho_fac_r(ijdim,kdim) !            rain
    real(8) :: rho_fac_i(ijdim,kdim) !            ice
    real(8) :: rho_fac_s(ijdim,kdim) !            snow 
    real(8) :: rho_fac_g(ijdim,kdim) !            graupel
    real(8) :: cva(ijdim,kdim)       ! 
    real(8) :: cpa(ijdim,kdim)       ! [Add] 09/08/18 T.Mitsui
    !
    real(8) :: drhogqv               ! d (rho*qv*gsgam2)
    real(8) :: drhogqc, drhognc      !        qc, nc
    real(8) :: drhogqr, drhognr      !        qr, nr
    real(8) :: drhogqi, drhogni      !        qi, ni
    real(8) :: drhogqs, drhogns      !        qs, ns
    real(8) :: drhogqg, drhogng      !        qg, ng

    ! production rate of nucleation
    real(8) :: PLCccn(ijdim,kdim), PNCccn(ijdim,kdim)
    real(8) :: PLIccn(ijdim,kdim), PNIccn(ijdim,kdim)
    ! production rate of freezing 
    real(8) :: PLChom(ijdim,kdim), PNChom(ijdim,kdim)
    real(8) :: PLChet(ijdim,kdim), PNChet(ijdim,kdim)
    real(8) :: PLRhet(ijdim,kdim), PNRhet(ijdim,kdim)
    ! production rate of melting
    real(8) :: PLImlt(ijdim,kdim), PNImlt(ijdim,kdim)
    real(8) :: PLSmlt(ijdim,kdim), PNSmlt(ijdim,kdim)
    real(8) :: PLGmlt(ijdim,kdim), PNGmlt(ijdim,kdim)
    ! production rate of vapor deposition
    real(8) :: PLRdep(ijdim,kdim), PNRdep(ijdim,kdim)
    real(8) :: PLIdep(ijdim,kdim), PNIdep(ijdim,kdim)
    real(8) :: PLSdep(ijdim,kdim), PNSdep(ijdim,kdim)
    real(8) :: PLGdep(ijdim,kdim), PNGdep(ijdim,kdim)
    real(8) :: PLCdep(ijdim,kdim), PNCdep(ijdim,kdim)

    ! production rate of warm collection process
    ! auto-conversion
    real(8) :: PLCaut(ijdim,kdim), PNCaut(ijdim,kdim)
    real(8) ::                     PNRaut(ijdim,kdim)
    ! accretion
    real(8) :: PLCacc(ijdim,kdim), PNCacc(ijdim,kdim)
    ! self-colletion, break-up
    real(8) :: PNRslc(ijdim,kdim), PNRbrk(ijdim,kdim)
    real(8) :: wrm_dqc, wrm_dnc
    real(8) :: wrm_dqr, wrm_dnr
    ! production rate of mixed-phase collection process
    ! PXXacYY2ZZ means XX collect YY produce ZZ
    real(8) :: PLIacLC2LI(ijdim,kdim)  , PNIacNC2NI(ijdim,kdim)   ! cloud-ice
    !
    real(8) :: PLSacLC2LS(ijdim,kdim),   PNSacNC2NS(ijdim,kdim)   ! cloud-snow(cloud change)
    !
    real(8) :: PLGacLC2LG(ijdim,kdim),   PNGacNC2NG(ijdim,kdim)   ! cloud-graupel
    real(8) :: PLRacLI2LG_I(ijdim,kdim), PNRacNI2NG_I(ijdim,kdim) ! rain-ice(ice change)
    real(8) :: PLRacLI2LG_R(ijdim,kdim), PNRacNI2NG_R(ijdim,kdim) ! rain-ice(rain change)
    real(8) :: PLRacLS2LG_S(ijdim,kdim), PNRacNS2NG_S(ijdim,kdim) ! rain-snow(snow change)
    real(8) :: PLRacLS2LG_R(ijdim,kdim), PNRacNS2NG_R(ijdim,kdim) ! rain-snow(rain change)
    real(8) :: PLRacLG2LG(ijdim,kdim),   PNRacNG2NG(ijdim,kdim)   ! rain-graupel(rain change)
    real(8) :: PLIacLI2LS(ijdim,kdim),   PNIacNI2NS(ijdim,kdim)   ! ice-ice
    real(8) :: PLIacLS2LS(ijdim,kdim),   PNIacNS2NS(ijdim,kdim)   ! ice-snow(ice change)
    real(8) :: PNSacNS2NS(ijdim,kdim)                             ! snow-snow
    real(8) :: PNGacNG2NG(ijdim,kdim)                             ! graupel-graupel 
    real(8) :: PLGacLS2LG(ijdim,kdim),   PNGacNS2NG(ijdim,kdim)   ! snow-graupel
    real(8) :: gc_dqc, gc_dnc
    real(8) :: sc_dqc, sc_dnc
    real(8) :: ic_dqc, ic_dnc
    real(8) :: rg_dqg, rg_dng
    real(8) :: rg_dqr, rg_dnr
    real(8) :: rs_dqr, rs_dnr, rs_dqs, rs_dns
    real(8) :: ri_dqr, ri_dnr
    real(8) :: ri_dqi, ri_dni
    real(8) :: ii_dqi, ii_dni
    real(8) :: is_dqi, is_dni, ss_dns
    real(8) :: gs_dqs, gs_dns, gg_dng 
    ! mixed-phase collection process total plus(clp_), total minus(clm_)
    real(8) :: clp_dqc, clp_dnc, clm_dqc, clm_dnc
    real(8) :: clp_dqr, clp_dnr, clm_dqr, clm_dnr
    real(8) :: clp_dqi, clp_dni, clm_dqi, clm_dni
    real(8) :: clp_dqs, clp_dns, clm_dqs, clm_dns
    real(8) :: clp_dqg, clp_dng, clm_dqg, clm_dng
    real(8) :: fac1, fac3, fac4, fac6, fac7, fac9, fac10
    ! production rate of partial conversion(ice, snow => graupel)
    real(8) :: PLIcon(ijdim,kdim), PNIcon(ijdim,kdim)
    real(8) :: PLScon(ijdim,kdim), PNScon(ijdim,kdim)
    real(8) :: pco_dqi, pco_dni
    real(8) :: pco_dqs, pco_dns
    real(8) :: pco_dqg, pco_dng
    ! production rate of enhanced melting due to 
    real(8) :: PLIacm(ijdim,kdim), PNIacm(ijdim,kdim) ! ice-cloud
    real(8) :: PLIarm(ijdim,kdim), PNIarm(ijdim,kdim) ! ice-rain
    real(8) :: PLSacm(ijdim,kdim), PNSacm(ijdim,kdim) ! snow-cloud
    real(8) :: PLSarm(ijdim,kdim), PNSarm(ijdim,kdim) ! snow-rain 
    real(8) :: PLGacm(ijdim,kdim), PNGacm(ijdim,kdim) ! graupel-cloud
    real(8) :: PLGarm(ijdim,kdim), PNGarm(ijdim,kdim) ! graupel-rain
    real(8) :: eml_dqc, eml_dnc
    real(8) :: eml_dqr, eml_dnr
    real(8) :: eml_dqi, eml_dni
    real(8) :: eml_dqs, eml_dns
    real(8) :: eml_dqg, eml_dng
    ! production rate of ice multiplication by splintering
    real(8) :: PLGspl(ijdim,kdim)
    real(8) :: PLSspl(ijdim,kdim)   
    real(8) :: PNIspl(ijdim,kdim)
    real(8) :: spl_dqi, spl_dni
    real(8) :: spl_dqg, spl_dqs

    real(8) :: rrhog(ijdim,kdim)
    !-----------------------------------------------
    ! work for explicit supersaturation modeling
    !-----------------------------------------------
    real(8) :: dTdt_equiv(ijdim,kdim)    !
    real(8) :: dlvsi(ijdim,kdim)
    !--------------------------------------------------
    !
    ! variables for output 
    !
    !--------------------------------------------------
    ! work for column production term                 
    real(8) :: sl_PLCdep(ijdim,1), sl_PNCdep(ijdim,1) ! 
    real(8) :: sl_PLCaut(ijdim,1), sl_PNCaut(ijdim,1) ! 
    real(8) :: sl_PLCacc(ijdim,1), sl_PNCacc(ijdim,1) ! 
    real(8) :: sl_PLRdep(ijdim,1), sl_PNRdep(ijdim,1) ! 

    real(8) :: ml_PNIcol(ijdim,kdim)
    real(8) :: ml_PLIcol(ijdim,kdim)
    real(8) :: ml_PNIfrz(ijdim,kdim) ! freezing
    real(8) :: ml_PNIspl(ijdim,kdim) ! ice-multiplication
    real(8) :: ml_PLCcnd(ijdim,kdim)  ! cnd+evp
    real(8) :: ml_PLRcnd(ijdim,kdim)  ! cnd+evp
    real(8) :: ml_PLIdep(ijdim,kdim)  ! dep+sbl
    real(8) :: ml_PLSdep(ijdim,kdim)  ! dep+sbl
    real(8) :: ml_PLGdep(ijdim,kdim)  ! dep+sbl
    !
    ! [Add] 10/08/03 T.Mitsui, heating rate 
    real(8) :: ml_dTacr(ijdim,kdim) ! accretion
    real(8) :: ml_dTcnd(ijdim,kdim)
    real(8) :: ml_dTdep(ijdim,kdim)
    real(8) :: ml_dTmlt(ijdim,kdim)
    real(8) :: ml_dTfrz(ijdim,kdim)
    real(8) :: ml_dTpc(ijdim,kdim)
    !--------------------------------------------------
    logical :: flag_history_in            ! 
    real(8) :: qke(ijdim,kdim)    !

    real(8), parameter :: eps=1.d-30
    real(8), parameter :: eps_qv=1.d-50    ! 
    real(8), parameter :: eps_rhoge=1.d-50 ! 
    real(8), parameter :: eps_rhog=1.d-50  ! 
    real(8) :: r_dt     
    real(8) :: wdt, r_wdt
    real(8) :: r_ntmax
    integer :: ntdiv


    real(8) :: rhoe(kdim,IA,JA)
    real(8) :: temp(kdim,IA,JA)
    real(8) :: pres(kdim,IA,JA)
    real(8) :: rhoq(kdim,IA,JA,nqmax)

    real(8) :: velw(kdim,IA,JA,nqmax)
    real(8) :: flux_rain (kdim,IA,JA)
    real(8) :: flux_snow (kdim,IA,JA)
    real(8) :: wflux_rain(kdim,IA,JA)
    real(8) :: wflux_snow(kdim,IA,JA)
    integer :: step

    integer :: k, i, j, iq, ij
    !---------------------------------------------------------------------------

    r_dt = 1.D0 / dt

    call TIME_rapstart('MP0 Setup')

    call negative_filter

    call TIME_rapend  ('MP0 Setup')
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

    do j = JS, JE
    do i = IS, IE
       ij = (j-JS)*IMAX+i-IS+1

       !OCL XFILL
       do k = kmin, kmax
          rhog  (ij,k) = var(k,i,j,I_DENS)
       enddo
       !OCL XFILL
       do k = kmin, kmax
          rhogw (ij,k) = var(k,i,j,I_MOMZ)
       enddo
       !OCL XFILL
       do k = kmin, kmax
          rhogvx(ij,k) = var(k,i,j,I_MOMX)
       enddo
       !OCL XFILL
       do k = kmin, kmax
          rhogvy(ij,k) = var(k,i,j,I_MOMY)
       enddo

       do k = kmin, kmax
          th    (ij,k) = var(k,i,j,I_RHOT) / var(k,i,j,I_DENS)
       enddo
    enddo
    enddo

    do ij = 1, ijdim
       rhog  (ij,1:kmin-1) = rhog  (ij,kmin)
       rhogw (ij,1:kmin-1) = rhogw (ij,kmin)
       rhogvx(ij,1:kmin-1) = rhogvx(ij,kmin)
       rhogvy(ij,1:kmin-1) = rhogvy(ij,kmin)
       th    (ij,1:kmin-1) = th    (ij,kmin)

       rhog  (ij,kmax+1:kdim) = rhog  (ij,kmax)
       rhogw (ij,kmax+1:kdim) = rhogw (ij,kmax)
       rhogvx(ij,kmax+1:kdim) = rhogvx(ij,kmax)
       rhogvy(ij,kmax+1:kdim) = rhogvy(ij,kmax)
       th    (ij,kmax+1:kdim) = th    (ij,kmax)
    enddo

    do iq = 1, nqmax
       do j = JS, JE
       do i = IS, IE
          ij = (j-JS)*IMAX+i-IS+1

          !OCL XFILL
          do k = kmin, kmax
             q(ij,k,iq) = var(k,i,j,5+iq)
          enddo
       enddo
       enddo
    enddo

    do iq = 1, nqmax
       do ij = 1, ijdim
          q(ij,     1:kmin-1,iq) = q(ij,kmin,iq)
          q(ij,kmax+1:kdim  ,iq) = q(ij,kmax,iq)
       enddo
    enddo

    do iq = 1, nqmax
    do k  = 1, kdim
    do ij = 1, ijdim
       rhogq(ij,k,iq) = rhog(ij,k) * q(ij,k,iq)
    enddo
    enddo
    enddo

    do k  = 1, kdim
    do ij = 1, ijdim
       rho  (ij,k) = rhog(ij,k)
       rrhog(ij,k) = 1.D0 / rhog(ij,k)
       w    (ij,k) = rhogw(ij,k) * rrhog(ij,k)
    enddo
    enddo

    call thrmdyn_qd( qd, q )
    call thrmdyn_cv( cva, q, qd )
    call thrmdyn_tempre2( tem, pre, rho, th, qd, q )

    do k  = 1, kdim
    do ij = 1, ijdim
       rhoge(ij,k) = rhog(ij,k) * tem(ij,k) * cva(ij,k)
    enddo
    enddo

    if( opt_debug_tem ) call debug_tem( 1, tem(:,:), rho(:,:), pre(:,:), q(:,:,I_QV) )

    call TIME_rapend  ('MPX ijkconvert')

    call TIME_rapstart('MP1 Nucleation')

    do k  = kmin, kmax
    do ij = 1,    ijdim
       rho_fac         = rho_0 / max(rho(ij,k),rho_min)
       rho_fac_c(ij,k) = rho_fac**gamma_v(I_QC)
       rho_fac_r(ij,k) = rho_fac**gamma_v(I_QR) 
       rho_fac_i(ij,k) = (pre(ij,k)/pre0_vt)**a_pre0_vt * (tem(ij,k)/tem0_vt)**a_tem0_vt
       rho_fac_s(ij,k) = rho_fac_i(ij,k)
       rho_fac_g(ij,k) = rho_fac_i(ij,k)
    enddo
    enddo
    !OCL XFILL
    do ij=1, ijdim
       rho_fac_c(ij,1:kmin-1)    = 1.D0
       rho_fac_r(ij,1:kmin-1)    = 1.D0
       rho_fac_i(ij,1:kmin-1)    = 1.D0
       rho_fac_s(ij,1:kmin-1)    = 1.D0
       rho_fac_g(ij,1:kmin-1)    = 1.D0

       rho_fac_c(ij,kmax+1:kdim) = rho_fac_c(ij,kmax)
       rho_fac_r(ij,kmax+1:kdim) = rho_fac_r(ij,kmax)
       rho_fac_i(ij,kmax+1:kdim) = rho_fac_i(ij,kmax)
       rho_fac_s(ij,kmax+1:kdim) = rho_fac_s(ij,kmax)
       rho_fac_g(ij,kmax+1:kdim) = rho_fac_g(ij,kmax)
    end do

    call thrmdyn_cp( cpa, q, qd )

    do k  = 1, kdim
    do ij = 1, ijdim
       wtem(ij,k) = max(tem(ij,k), tem_min) 
    enddo
    enddo

    !OCL XFILL
    do k  = 1, kdim
    do ij = 1, ijdim
       qke       (ij,k) = 0.D0 ! 2*TKE 
       dTdt_equiv(ij,k) = 0.D0
    enddo
    enddo

    do k  = 1, kdim
    do ij = 1, ijdim
       lv(ij,k) = rho(ij,k)*q(ij,k,I_QV)
       ni(ij,k) = max( 0.D0, rho(ij,k)*q(ij,k,I_NI) )
       nc(ij,k) = max( 0.D0, rho(ij,k)*q(ij,k,I_NC) )
    enddo
    enddo

    call nucleation(     &
         ijdim, kdim,    & ! in
         kmin, kmax,     & ! in
         l_region,       & ! in 
         z, w,           & ! in 
         rho, wtem, pre, & ! in
         lv, nc, ni,     & ! in
         PNCccn, PLCccn, & ! out
         PNIccn, PLIccn, & ! out
         cva, cpa,       & ! in 
         dTdt_equiv,     & ! in 
         qke,            & ! in 
         flag_history_in,& ! in 
         ct, dt         ) ! in 

    do k  = kmin, kmax
    do ij = 1,    ijdim
       ! nucleation
       drhogqc = dt * PLCccn(ij,k)
       drhognc = dt * PNCccn(ij,k)
       drhogqi = dt * PLIccn(ij,k)
       drhogni = dt * PNIccn(ij,k)
       drhogqv = max( -rhogq(ij,k,I_QV), -drhogqc-drhogqi )
       fac1    = drhogqv / min( -drhogqc-drhogqi, -eps ) ! limiting coefficient

       rhogq(ij,k,I_QV) = rhogq(ij,k,I_QV) + drhogqv
       rhogq(ij,k,I_QC) = max(0.D0, rhogq(ij,k,I_QC) + drhogqc*fac1)
       rhogq(ij,k,I_QI) = max(0.D0, rhogq(ij,k,I_QI) + drhogqi*fac1)
       rhogq(ij,k,I_NC) = max(0.D0, rhogq(ij,k,I_NC) + drhognc)
       rhogq(ij,k,I_NI) = max(0.D0, rhogq(ij,k,I_NI) + drhogni)

       ! cloud number concentration filter
       rhogq(ij,k,I_NC) = min( rhogq(ij,k,I_NC), nc_uplim(ij,1) )
    enddo
    enddo

    do k  = kmin, kmax
    do ij = 1,    ijdim
       rhoge(ij,k) = rhoge(ij,k) - LHV * drhogqv + LHF * drhogqi*fac1
    enddo
    enddo

    do k  = kmin, kmax
    do ij = 1,    ijdim
       q(ij,k,I_QV) = rhogq(ij,k,I_QV) * rrhog(ij,k)
       q(ij,k,I_QC) = rhogq(ij,k,I_QC) * rrhog(ij,k)
       q(ij,k,I_QI) = rhogq(ij,k,I_QI) * rrhog(ij,k)
       q(ij,k,I_NC) = rhogq(ij,k,I_NC) * rrhog(ij,k)
       q(ij,k,I_NI) = rhogq(ij,k,I_NI) * rrhog(ij,k)
    enddo
    enddo

    !--- update mass concentration  
    call thrmdyn_qd( qd, q )
    call thrmdyn_cv( cva, q, qd )

    do k  = kmin, kmax
    do ij = 1,    ijdim
       tem(ij,k) = rhoge(ij,k) / ( rhog(ij,k) * cva(ij,k) )
       pre(ij,k) = rho(ij,k)*( qd(ij,k)*CNST_RAIR+q(ij,k,I_QV)*CNST_RVAP )*tem(ij,k)
    enddo
    enddo

    if( opt_debug )     call debugreport_nucleation
    if( opt_debug_tem ) call debug_tem( 2, tem(:,:), rho(:,:), pre(:,:), q(:,:,I_QV) )

    call TIME_rapend  ('MP1 Nucleation')
    !----------------------------------------------------------------------------
    !
    ! 2.Phase change: Freezing, Melting, Vapor deposition
    ! 
    !----------------------------------------------------------------------------
    call TIME_rapstart('MP2 Phase change')

    ! parameter setting
    wdt=dt
    r_wdt=1.d0/wdt

       if(  ntdiv     == ntmax_phase_change  )then
          flag_history_in=.true.
       end if
       !
       do k=1, kdim
          do ij=1, ijdim
             lr(ij,k)     = rhogq(ij,k,I_QR)
             nr(ij,k)     = rhogq(ij,k,I_NR)
             xr(ij,k)     = max(xr_min,  min(xr_max, lr(ij,k)/(nr(ij,k)+nr_min) ))
             dr_xa(ij,k)  = a_m(I_QR)*xr(ij,k)**b_m(I_QR)
             vt_xa(ij,k,I_QR,1) = alpha_v(I_QR,1)*(xr(ij,k)**beta_v(I_QR,1))*rho_fac_r(ij,k)
             vt_xa(ij,k,I_QR,2) = vt_xa(ij,k,I_QR,1)
          end do
       end do
       !
       do k=1, kdim
          do ij=1, ijdim
             !! Following values shoud be already filtered to be non-zero before sbroutine was called.
             ! Mass concentration [kg/m3]
             lv(ij,k)     = rhogq(ij,k,I_QV)
             !
             lc(ij,k)     = rhogq(ij,k,I_QC)
             li(ij,k)     = rhogq(ij,k,I_QI)
             ls(ij,k)     = rhogq(ij,k,I_QS)
             lg(ij,k)     = rhogq(ij,k,I_QG)
             ! Number concentration[/m3] (should be filtered to avoid zero division.)
             nc(ij,k)     = rhogq(ij,k,I_NC)
             ni(ij,k)     = rhogq(ij,k,I_NI)
             ns(ij,k)     = rhogq(ij,k,I_NS)
             ng(ij,k)     = rhogq(ij,k,I_NG)
             ! Mass of mean particle [kg]  
             ! SB06(94)
             !
             xc(ij,k)     = min(xc_max, max(xc_min, lc(ij,k)/(nc(ij,k)+nc_min) ))
             xi(ij,k)     = min(xi_max, max(xi_min, li(ij,k)/(ni(ij,k)+ni_min) ))
             xs(ij,k)     = min(xs_max, max(xs_min, ls(ij,k)/(ns(ij,k)+ns_min) ))
             xg(ij,k)     = min(xg_max, max(xg_min, lg(ij,k)/(ng(ij,k)+ng_min) ))
             ! diamter of average mass
             ! SB06(32)
             dc_xa(ij,k)  = a_m(I_QC)*xc(ij,k)**b_m(I_QC)
             di_xa(ij,k)  = a_m(I_QI)*xi(ij,k)**b_m(I_QI)
             ds_xa(ij,k)  = a_m(I_QS)*xs(ij,k)**b_m(I_QS)
             dg_xa(ij,k)  = a_m(I_QG)*xg(ij,k)**b_m(I_QG)
             ! terminal velocity of average mass
             vt_xa(ij,k,I_QC,1) = alpha_v(I_QC,1)*(xc(ij,k)**beta_v(I_QC,1))*rho_fac_c(ij,k)
             vt_xa(ij,k,I_QI,1) = alpha_v(I_QI,1)*(xi(ij,k)**beta_v(I_QI,1))*rho_fac_i(ij,k)
             vt_xa(ij,k,I_QS,1) = alpha_v(I_QS,1)*(xs(ij,k)**beta_v(I_QS,1))*rho_fac_s(ij,k)
             vt_xa(ij,k,I_QG,1) = alpha_v(I_QG,1)*(xg(ij,k)**beta_v(I_QG,1))*rho_fac_g(ij,k)
             vt_xa(ij,k,I_QC,2) = alpha_v(I_QC,2)*(xc(ij,k)**beta_v(I_QC,2))*rho_fac_c(ij,k)
             vt_xa(ij,k,I_QI,2) = alpha_v(I_QI,2)*(xi(ij,k)**beta_v(I_QI,2))*rho_fac_i(ij,k)
             vt_xa(ij,k,I_QS,2) = alpha_v(I_QS,2)*(xs(ij,k)**beta_v(I_QS,2))*rho_fac_s(ij,k)
             vt_xa(ij,k,I_QG,2) = alpha_v(I_QG,2)*(xg(ij,k)**beta_v(I_QG,2))*rho_fac_g(ij,k)
          end do
       end do
       !
       wtem(:,:) = max(tem(:,:), tem_min) 
       call moist_psat_water( wtem, esw ) 
       call moist_psat_ice( wtem, esi )
       !
       if( flag_history_in )then
          dlvsi(:,:) = lv(:,:) - esi(:,:)/(CNST_RVAP*tem(:,:))
       end if

       call thrmdyn_qd( qd, q )

       call freezing_water( &
            ijdim, kdim,    & ! in
            kmin, kmax,     & ! in
            wdt,            & ! in 
            PLChom, PNChom, & ! out
            PLChet, PNChet, & ! out
            PLRhet, PNRhet, & ! out
            lc, lr, nc, nr, xc, xr, tem ) ! in 
       !
       call dep_vapor_melt_ice( &
            ijdim, kdim,      & ! in
            kmin, kmax,       & ! in
            PLCdep,           & ! out
            PLRdep, PNRdep,   & ! out
            PLIdep, PNIdep,   & ! out
            PLSdep, PNSdep,   & ! out
            PLGdep, PNGdep,   & ! out
            PLImlt, PNImlt,   & ! out
            PLSmlt, PNSmlt,   & ! out
            PLGmlt, PNGmlt,   & ! out
            rho, wtem, pre, lv,& ! in 
            qd,               & ! in
            esw, esi,         & ! in
            nc,               & ! in
            nr, ni, ns, ng,   & ! in
            li, ls, lg,       & ! in
            xc,               & ! in
            xr, xi, xs, xg,   & ! in
            vt_xa,            & ! in
            dc_xa,            & ! in
            dr_xa, di_xa,     & ! in
            ds_xa, dg_xa      ) ! in
       !
       ! update subroutine
       !    
       call update_by_phase_change(   &
            ntdiv, ntmax_phase_change,& ! in 
            ijdim, kdim, kmin, kmax,  & ! in
            nqmax,                    & ! in
            wdt,                      & ! in
            gsgam2,                   & ! in
            z,                        & ! in
            dz,                       & ! in
            w,                        & ! in
            dTdt_equiv,               & ! in 
            rhog,                     & ! in
            rhoge,                    & ! inout
            rhogq(:,:,:), q(:,:,:),   & ! inout
            tem, pre, rho,            & ! inout
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
            PLCdep,         & ! inout
            PLRdep, PNRdep, & ! inout
            PLIdep, PNIdep, & ! inout
            PLSdep, PNSdep, & ! inout
            PLGdep, PNGdep, & ! inout
            ! melting term
            PLImlt, PNImlt, & ! in
            PLSmlt, PNSmlt, & ! in
            PLGmlt, PNGmlt, & ! in
            flag_history_in,& ! in
            sl_PLCdep,      &      ! out
            sl_PLRdep, sl_PNRdep ) ! out

       if( opt_debug )     call debugreport_phasechange
       if( opt_debug_tem ) call debug_tem( 3, tem(:,:), rho(:,:), pre(:,:), q(:,:,I_QV) )

    call TIME_rapend  ('MP2 Phase change')
    !----------------------------------------------------------------------------
    !
    ! 3.Collection process
    ! 
    !----------------------------------------------------------------------------
    call TIME_rapstart('MP3 Collection')

    wdt = dt
    flag_history_in=.true.

    ! parameter setting
    do k  = 1, kdim
    do ij = 1, ijdim
       ! Mass concentration [kg/m3]
       lv(ij,k) = rhogq(ij,k,I_QV) 
       lc(ij,k) = rhogq(ij,k,I_QC) 
       lr(ij,k) = rhogq(ij,k,I_QR) 
       li(ij,k) = rhogq(ij,k,I_QI) 
       ls(ij,k) = rhogq(ij,k,I_QS) 
       lg(ij,k) = rhogq(ij,k,I_QG) 
    enddo
    enddo
    do k  = 1, kdim
    do ij = 1, ijdim
       ! Number concentration[/m3]
       nc(ij,k) = rhogq(ij,k,I_NC) 
       nr(ij,k) = rhogq(ij,k,I_NR) 
       ni(ij,k) = rhogq(ij,k,I_NI) 
       ns(ij,k) = rhogq(ij,k,I_NS) 
       ng(ij,k) = rhogq(ij,k,I_NG) 
    enddo
    enddo
    do k  = 1, kdim
    do ij = 1, ijdim
       ! Mass of mean particle [kg]
       xc(ij,k) = min(xc_max, max(xc_min, lc(ij,k)/(nc(ij,k)+nc_min) ) )
       xr(ij,k) = min(xr_max, max(xr_min, lr(ij,k)/(nr(ij,k)+nr_min) ) )
       xi(ij,k) = min(xi_max, max(xi_min, li(ij,k)/(ni(ij,k)+ni_min) ) )
       xs(ij,k) = min(xs_max, max(xs_min, ls(ij,k)/(ns(ij,k)+ns_min) ) )
       xg(ij,k) = min(xg_max, max(xg_min, lg(ij,k)/(ng(ij,k)+ng_min) ) )
    enddo
    enddo
    do k  = 1, kdim
    do ij = 1, ijdim
       ! effective cross section is assume as area equivalent circle 
       dc_xa(ij,k)  = 2.d0*a_rea(I_QC)*xc(ij,k)**b_rea(I_QC)
       dr_xa(ij,k)  = 2.d0*a_rea(I_QR)*xr(ij,k)**b_rea(I_QR)
       di_xa(ij,k)  = 2.d0*a_rea(I_QI)*xi(ij,k)**b_rea(I_QI)
       ds_xa(ij,k)  = 2.d0*a_rea(I_QS)*xs(ij,k)**b_rea(I_QS)
       dg_xa(ij,k)  = 2.d0*a_rea(I_QG)*xg(ij,k)**b_rea(I_QG)
    enddo
    enddo
    do k  = 1, kdim
    do ij = 1, ijdim
       ! terminal velocity of average mass
       ! SB06(33)
       vt_xa(ij,k,I_QC,2) = alpha_v(I_QC,2)*(xc(ij,k)**beta_v(I_QC,2))*rho_fac_c(ij,k)
       vt_xa(ij,k,I_QR,2) = alpha_v(I_QR,2)*(xr(ij,k)**beta_v(I_QR,2))*rho_fac_r(ij,k)
       vt_xa(ij,k,I_QI,2) = alpha_v(I_QI,2)*(xi(ij,k)**beta_v(I_QI,2))*rho_fac_i(ij,k)
       vt_xa(ij,k,I_QS,2) = alpha_v(I_QS,2)*(xs(ij,k)**beta_v(I_QS,2))*rho_fac_s(ij,k)
       vt_xa(ij,k,I_QG,2) = alpha_v(I_QG,2)*(xg(ij,k)**beta_v(I_QG,2))*rho_fac_g(ij,k)
    enddo
    enddo

       call aut_acc_slc_brk(  &
            ijdim, kdim,      &
            kmin, kmax,       &
            PLCaut, PNCaut,   &
            PNRaut,           &
            PLCacc, PNCacc,   &
            PNRslc, PNRbrk,   &
            lc, lr, nc, nr,xc,&
            dr_xa,            &
            rho, tem          )
       !
       call mixed_phase_collection(         &
            ijdim, kdim,                    & ! in
            kmin, kmax,                     & ! in
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
            PNGacNG2NG,                     & ! out
            PLGacLS2LG,   PNGacNS2NG,       & ! out
            ! partial conversion (ice, snow => graupel)
            PLIcon, PNIcon, PLScon, PNScon, & ! out
            ! enhanced melting (latent heat effect )
            PLIacm, PNIacm, PLIarm, PNIarm, & ! out
            PLSacm, PNSacm, PLSarm, PNSarm, & ! out
            PLGacm, PNGacm, PLGarm, PNGarm, & ! out
            tem,                            & ! in 
            lv,                             & ! in
            lc, lr, li, ls, lg,             & ! in
            nc, nr, ni, ns, ng,             & ! in
            xc, xr, xi, xs, xg,             & ! in
            dc_xa, dr_xa, di_xa, ds_xa, dg_xa,& ! in
            vt_xa(:,:,I_QC,2), vt_xa(:,:,I_QR,2), &
            vt_xa(:,:,I_QI,2), vt_xa(:,:,I_QS,2), vt_xa(:,:,I_QG,2), &
            rho(:,:),                       & ! in 
            flag_history_in                 ) ! in 
       !
       call ice_multiplication( &
            ijdim, kdim,        & ! in
            kmin, kmax,         & ! in
            PLGspl, PLSspl,     & ! out
            PNIspl,             & ! out
            PLIacLC2LI,         & ! in
            PLSacLC2LS,         & ! in
            PLGacLC2LG,         & ! in
            tem, lc, nc, xi, xc ) ! in
       !
       ! update
       ! rhogq = l*gsgam
       !
       do k=kmin, kmax
          do ij=1, ijdim
             ! 
             ! warm collection process 
             wrm_dqc = max( wdt*( PLCaut(ij,k)+PLCacc(ij,k) ), -lc(ij,k)  )
             wrm_dnc = max( wdt*( PNCaut(ij,k)+PNCacc(ij,k) ), -nc(ij,k)  )
             wrm_dnr = max( wdt*( PNRaut(ij,k)+PNRslc(ij,k)+PNRbrk(ij,k) ), -nr(ij,k) )
             wrm_dqr = -wrm_dqc
             ! mixed phase collection
             ! Pxxacyy2zz xx and yy decrease and zz increase .
             !
             ! At first fixer is applied to decreasing particles.
             ! order of fixer: graupel-cloud, snow-cloud, ice-cloud, graupel-rain, snow-rain, ice-rain, 
             !                 snow-ice,  ice-ice, graupel-snow, snow-snow
             ! cloud mass decrease 
             gc_dqc  = max( wdt*PLGacLC2LG(ij,k)  , min(0.D0, -lc(ij,k)-wrm_dqc               )) ! => dqg
             sc_dqc  = max( wdt*PLSacLC2LS(ij,k)  , min(0.D0, -lc(ij,k)-wrm_dqc-gc_dqc        )) ! => dqs
             ic_dqc  = max( wdt*PLIacLC2LI(ij,k)  , min(0.D0, -lc(ij,k)-wrm_dqc-gc_dqc-sc_dqc )) ! => dqi
             ! cloud num. decrease 
             gc_dnc  = max( wdt*PNGacNC2NG(ij,k)  , min(0.D0, -nc(ij,k)-wrm_dnc               )) ! => dnc:minus
             sc_dnc  = max( wdt*PNSacNC2NS(ij,k)  , min(0.D0, -nc(ij,k)-wrm_dnc-gc_dnc        )) ! => dnc:minus
             ic_dnc  = max( wdt*PNIacNC2NI(ij,k)  , min(0.D0, -nc(ij,k)-wrm_dnc-gc_dnc-sc_dnc )) ! => dnc:minus
             ! rain mass decrease ( tem < 273.15K)
             if( tem(ij,k) <= CNST_TEM00 )then
                rg_dqr  = max( wdt*PLRacLG2LG  (ij,k), min(0.D0, -lr(ij,k)-wrm_dqr               )) ! => dqg
                rg_dqg  = 0.D0
                rs_dqr  = max( wdt*PLRacLS2LG_R(ij,k), min(0.D0, -lr(ij,k)-wrm_dqr-rg_dqr        )) ! => dqg
                ri_dqr  = max( wdt*PLRacLI2LG_R(ij,k), min(0.D0, -lr(ij,k)-wrm_dqr-rg_dqr-rs_dqr )) ! => dqg
                ! rain num. decrease 
                rg_dnr  = max( wdt*PNRacNG2NG  (ij,k), min(0.D0, -nr(ij,k)-wrm_dnr               )) ! => dnr:minus,dng:plus
                rg_dng  = 0.D0
                rs_dnr  = max( wdt*PNRacNS2NG_R(ij,k), min(0.D0, -nr(ij,k)-wrm_dnr-rg_dnr        )) ! => dnr:minus,dng:plus
                ri_dnr  = max( wdt*PNRacNI2NG_R(ij,k), min(0.D0, -nr(ij,k)-wrm_dnr-rg_dnr-rs_dnr )) ! => dnr:minus,dng:plus
             else
                rg_dqg  = max( wdt*PLRacLG2LG  (ij,k), min(0.D0, -lg(ij,k)                       )) ! => dqg
                rg_dqr  = 0.D0 ! r+g -> r
                rs_dqr  = 0.D0 ! r+s -> r
                ri_dqr  = 0.D0 ! r+i -> r
                ! rain num. decrease 
                rg_dng  = max( wdt*PNRacNG2NG  (ij,k), min(0.D0, -ng(ij,k)                       )) ! => dnr:minus,dng:plus
                rg_dnr  = 0.D0 ! r+g -> r
                rs_dnr  = 0.D0 ! r+s -> r
                ri_dnr  = 0.D0 ! r+i -> r
             end if
             ! ice mass decrease 
             fac1    = (ri_dqr-eps)/ (wdt*PLRacLI2LG_R(ij,k)-eps) ! suppress factor by filter of rain
             ri_dqi  = max( wdt*PLRacLI2LG_I(ij,k)*fac1, min(0.D0, -li(ij,k)+ic_dqc               )) ! => dqg
             ii_dqi  = max( wdt*PLIacLI2LS(ij,k)       , min(0.D0, -li(ij,k)+ic_dqc-ri_dqi        )) ! => dqs
             is_dqi  = max( wdt*PLIacLS2LS(ij,k)       , min(0.D0, -li(ij,k)+ic_dqc-ri_dqi-ii_dqi )) ! => dqs
             ! ice num. decrease 
             fac4    = (ri_dnr-eps)/ (wdt*PNRacNI2NG_R(ij,k)-eps) ! suppress factor by filter of rain
             ri_dni  = max( wdt*PNRacNI2NG_I(ij,k)*fac4, min(0.D0, -ni(ij,k)               )) ! => dni:minus
             ii_dni  = max( wdt*PNIacNI2NS(ij,k)       , min(0.D0, -ni(ij,k)-ri_dni        )) ! => dni:minus,dns:plus(*0.5)
             is_dni  = max( wdt*PNIacNS2NS(ij,k)       , min(0.D0, -ni(ij,k)-ri_dni-ii_dni )) ! => dni:minus,dns:plus
             ! snow mass decrease 
             fac3    = (rs_dqr-eps)/(wdt*PLRacLS2LG_R(ij,k)-eps) ! suppress factor by filter of rain
             rs_dqs  = max( wdt*PLRacLS2LG_S(ij,k)*fac3, min(0.D0, -ls(ij,k)+sc_dqc+ii_dqi+is_dqi        )) ! => dqg
             gs_dqs  = max( wdt*PLGacLS2LG(ij,k)       , min(0.D0, -ls(ij,k)+sc_dqc+ii_dqi+is_dqi-rs_dqs )) ! => dqg
             ! snow num. decrease 
             fac6    = (rs_dnr-eps)/(wdt*PNRacNS2NG_R(ij,k)-eps) ! suppress factor by filter of rain
             fac7    = (is_dni-eps)/(wdt*PNIacNS2NS  (ij,k)-eps) ! suppress factor by filter of ice
             rs_dns  = max( wdt*PNRacNS2NG_S(ij,k)*fac6, min(0.D0, -ns(ij,k)+0.5d0*ii_dni+is_dni       )) ! => dns:minus
             gs_dns  = max( wdt*PNGacNS2NG(ij,k)       , min(0.D0, -ns(ij,k)+0.5d0*ii_dni+is_dni-rs_dns )) ! => dns:minus
             ss_dns  = max( wdt*PNSacNS2NS(ij,k)       , min(0.D0, -ns(ij,k)+0.5d0*ii_dni+is_dni-rs_dns-gs_dns ))
             !
             gg_dng  = max( wdt*PNGacNG2NG(ij,k)       , min(0.D0, -ng(ij,k) ))
             !          
             ! total plus in mixed phase collection(clp_)
             ! mass
             if( tem(ij,k) <= CNST_TEM00 )then
                clp_dqc =  0.D0
                clp_dqr =  0.D0
                clp_dqi = -ic_dqc
                clp_dqs = -sc_dqc-ii_dqi-is_dqi
                clp_dqg = -gc_dqc-rg_dqr-rs_dqr-rs_dqs-ri_dqr-ri_dqi-gs_dqs
                ! num.( number only increase when a+b=>c,  dnc=-dna)
                clp_dnc = 0.D0
                clp_dnr = 0.D0
                clp_dni = 0.D0
                clp_dns = -ii_dni*0.5d0
                clp_dng =       -rs_dnr-ri_dnr
                ! total minus in mixed phase collection(clm_)
                ! mass
                clm_dqc = gc_dqc+sc_dqc+ic_dqc
                clm_dqr = rg_dqr+rs_dqr+ri_dqr
                clm_dqi = ri_dqi+ii_dqi+is_dqi
                clm_dqs = rs_dqs+gs_dqs
                clm_dqg = 0.D0
                ! num.
                clm_dnc = gc_dnc+sc_dnc+ic_dnc
                clm_dnr = rg_dnr+rs_dnr+ri_dnr
                clm_dni = ri_dni+ii_dni+is_dni
                clm_dns = rs_dns+ss_dns+gs_dns
!!$             clm_dng = 0.D0
                clm_dng = gg_dng ! [Mod] 11/08/30 T.Mitsui
             else
                clp_dqc =  0.D0
                clp_dqr = -rg_dqg-rs_dqs-ri_dqi
                clp_dqi = -ic_dqc
                clp_dqs = -sc_dqc-ii_dqi-is_dqi
                clp_dqg = -gc_dqc-gs_dqs
                ! num.( number only increase when a+b=>c,  dnc=-dna)
                clp_dnc = 0.D0
                clp_dnr = 0.D0
                clp_dni = 0.D0
                clp_dns = -ii_dni*0.5d0
                clp_dng = 0.D0
                ! total minus in mixed phase collection(clm_)
                ! mass
                clm_dqc = gc_dqc+sc_dqc+ic_dqc
                clm_dqr = 0.D0
                clm_dqi = ri_dqi+ii_dqi+is_dqi
                clm_dqs = rs_dqs+gs_dqs
                clm_dqg = rg_dqg
                ! num.
                clm_dnc = gc_dnc+sc_dnc+ic_dnc
                clm_dnr = 0.D0
                clm_dni = ri_dni+ii_dni+is_dni
                clm_dns = rs_dns+ss_dns+gs_dns
!!$             clm_dng = rg_dng
                clm_dng = rg_dng+gg_dng ! [Mod] 11/08/30 T.Mitsui
             end if
             !
             ! partial conversion
             ! 08/05/08 [Mod] T.Mitsui
             pco_dqi = max( wdt*PLIcon(ij,k), -clp_dqi )
             pco_dqs = max( wdt*PLScon(ij,k), -clp_dqs )
             pco_dqg = -pco_dqi-pco_dqs
             ! 08/05/08 [Mod] T.Mitsui
             pco_dni = max( wdt*PNIcon(ij,k), -clp_dni )
             pco_dns = max( wdt*PNScon(ij,k), -clp_dns )
             pco_dng = -pco_dni-pco_dns
             ! enhanced melting ( always negative value )
             ! ice-cloud melting produces cloud, others produce rain
             eml_dqi =  max( wdt*PLIacm(ij,k), min(0.D0, -li(ij,k)-(clp_dqi+clm_dqi)-pco_dqi ))
             eml_dqs =  max( wdt*PLSacm(ij,k), min(0.D0, -ls(ij,k)-(clp_dqs+clm_dqs)-pco_dqs ))
             eml_dqg =  max( wdt*(PLGacm(ij,k)+PLGarm(ij,k)+PLSarm(ij,k)+PLIarm(ij,k)), &
                  min(0.D0, -lg(ij,k)-(clp_dqg+clm_dqg)-pco_dqg ))
             eml_dqc = -eml_dqi
             eml_dqr = -eml_dqs-eml_dqg
             !
             eml_dni =  max( wdt*PNIacm(ij,k), min(0.D0, -ni(ij,k)-(clp_dni+clm_dni)-pco_dni ))
             eml_dns =  max( wdt*PNSacm(ij,k), min(0.D0, -ns(ij,k)-(clp_dns+clm_dns)-pco_dns ))
             eml_dng =  max( wdt*(PNGacm(ij,k)+PNGarm(ij,k)+PNSarm(ij,k)+PNIarm(ij,k)), &
                  min(0.D0, -ng(ij,k)-(clp_dng+clm_dng)-pco_dng ))
             eml_dnc = -eml_dni
             eml_dnr = -eml_dns-eml_dng
             !
             ! ice multiplication
             spl_dqg = max( wdt*PLGspl(ij,k), min(0.D0, -lg(ij,k)-(clp_dqg+clm_dqg)-pco_dqg-eml_dqg ))
             spl_dqs = max( wdt*PLSspl(ij,k), min(0.D0, -ls(ij,k)-(clp_dqs+clm_dqs)-pco_dqs-eml_dqs ))
             spl_dqi = -spl_dqg-spl_dqs
             fac9    = (spl_dqg-eps)/(wdt*PLGspl(ij,k)-eps)
             fac10   = (spl_dqs-eps)/(wdt*PLSspl(ij,k)-eps)
             spl_dni = wdt*PNIspl(ij,k)*fac9*fac10
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
             rhogq(ij,k,I_QC) = max(0.D0, rhogq(ij,k,I_QC) + drhogqc )
             rhogq(ij,k,I_NC) = max(0.D0, rhogq(ij,k,I_NC) + drhognc )
             rhogq(ij,k,I_QR) = max(0.D0, rhogq(ij,k,I_QR) + drhogqr )
             rhogq(ij,k,I_NR) = max(0.D0, rhogq(ij,k,I_NR) + drhognr )
             rhogq(ij,k,I_QI) = max(0.D0, rhogq(ij,k,I_QI) + drhogqi )
             rhogq(ij,k,I_NI) = max(0.D0, rhogq(ij,k,I_NI) + drhogni )        
             rhogq(ij,k,I_QS) = max(0.D0, rhogq(ij,k,I_QS) + drhogqs )
             rhogq(ij,k,I_NS) = max(0.D0, rhogq(ij,k,I_NS) + drhogns )        
             rhogq(ij,k,I_QG) = max(0.D0, rhogq(ij,k,I_QG) + drhogqg )
             rhogq(ij,k,I_NG) = max(0.D0, rhogq(ij,k,I_NG) + drhogng )
          end do
       end do

    do k  = kmin, kmax
    do ij = 1,    ijdim
       rhoge(ij,k) = rhoge(ij,k) + LHF * ( drhogqi + drhogqs + drhogqg )
    enddo
    enddo

    !--- update mixing ratio
    do iq = 1,    nqmax
    do k  = kmin, kmax
    do ij = 1,    ijdim
       q(ij,k,iq) = rhogq(ij,k,iq) * rrhog(ij,k)
    enddo
    enddo
    enddo

    call thrmdyn_qd( qd, q )
    call thrmdyn_cv( cva, q, qd )

    do k  = kmin, kmax
    do ij = 1,    ijdim
       tem(ij,k) = rhoge(ij,k) / ( rhog(ij,k) * cva(ij,k) )
       pre(ij,k) = rho(ij,k)*( qd(ij,k)*CNST_RAIR+q(ij,k,I_QV)*CNST_RVAP )*tem(ij,k)
    enddo
    enddo

    if( opt_debug )     call debugreport_collection
    if( opt_debug_tem ) call debug_tem( 4, tem(:,:), rho(:,:), pre(:,:), q(:,:,I_QV) )

    call TIME_rapend  ('MP3 Collection')

    call TIME_rapstart('MPX ijkconvert')

    do k  = kmin, kmax
    do ij = 1,    ijdim
       th(ij,k) = tem(ij,k) * ( PRE00 / pre(ij,k) )**RovCP
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       ij = (j-JS)*IMAX+i-IS+1

       !OCL XFILL
       do k = kmin, kmax
          var(k,i,j,I_DENS) = rhog  (ij,k)
       enddo
       !OCL XFILL
       do k = kmin, kmax
          var(k,i,j,I_MOMZ) = rhogw (ij,k)
       enddo
       !OCL XFILL
       do k = kmin, kmax
          var(k,i,j,I_MOMX) = rhogvx(ij,k)
       enddo
       !OCL XFILL
       do k = kmin, kmax
          var(k,i,j,I_MOMY) = rhogvy(ij,k)
       enddo

       do k = kmin, kmax
          var(k,i,j,I_RHOT) = th    (ij,k) * rhog(ij,k)
       enddo

       !OCL XFILL
       do k = kmin, kmax
          rhoe(k,i,j) = rhoge(ij,k)
       enddo
       !OCL XFILL
       do k = kmin, kmax
          temp(k,i,j) = tem  (ij,k)
       enddo
       !OCL XFILL
       do k = kmin, kmax
          pres(k,i,j) = pre  (ij,k)
       enddo
    enddo
    enddo

    !OCL SERIAL
    do iq = 1, nqmax
       !OCL PARALLEL
       do j  = JS, JE
       do i  = IS, IE
          ij = (j-JS)*IMAX+i-IS+1

          !OCL XFILL
          do k = kmin, kmax
             rhoq(k,i,j,iq) = rhogq(ij,k,iq)
          enddo
       enddo
       enddo

       !OCL PARALLEL
       do j = JS, JE
       do i = IS, IE
       do k = kmin, kmax
          var(k,i,j,5+iq) = rhoq(k,i,j,iq) / var(k,i,j,I_DENS)
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

    !OCL XFILL
    do j = JS, JE
    do i = IS, IE
    do k = kmin, kmax
       flux_rain(k,i,j) = 0.D0
       flux_snow(k,i,j) = 0.D0
    enddo
    enddo
    enddo

    do step = 1, MP_NSTEP_SEDIMENTATION 

       call terminal_velocity( velw(:,:,:,:), &
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
       do k = kmin, kmax
          flux_rain(k,i,j) = flux_rain(k,i,j) + wflux_rain(k,i,j) * MP_RNSTEP_SEDIMENTATION
          flux_snow(k,i,j) = flux_snow(k,i,j) + wflux_snow(k,i,j) * MP_RNSTEP_SEDIMENTATION
       enddo
       enddo
       enddo

       if( opt_debug ) call debugreport_sedimentation

    enddo

    call TIME_rapend  ('MP5 Sedimentation')


    !----------------------------------------------------------------------------
    ! 6.Filter for rounding error(negative value) and artificial number
    !   ( We assume artificial filter as evaporation of the smallest particles )
    !----------------------------------------------------------------------------
    call TIME_rapstart('MP6 Filter')

    call negative_filter

    call TIME_rapend  ('MP6 Filter')

    return
    !---------------------------------------------------------------------------
  contains
    !---------------------------------------------------------------------------
    subroutine debugreport_nucleation
      implicit none

      write(IO_FID_LOG,*) "*** Nucleation/mod_mp_ndw6"

      write(IO_FID_LOG,*) "I_QV: max",maxval(rhogq(:,:,I_QV))," min:",minval(rhogq(:,:,I_QV))
      write(IO_FID_LOG,*) "I_QC: max",maxval(rhogq(:,:,I_QC))," min:",minval(rhogq(:,:,I_QC))
      write(IO_FID_LOG,*) "I_QR: max",maxval(rhogq(:,:,I_QR))," min:",minval(rhogq(:,:,I_QR))
      write(IO_FID_LOG,*) "I_QI: max",maxval(rhogq(:,:,I_QI))," min:",minval(rhogq(:,:,I_QI))
      write(IO_FID_LOG,*) "I_QS: max",maxval(rhogq(:,:,I_QS))," min:",minval(rhogq(:,:,I_QS))
      write(IO_FID_LOG,*) "I_QG: max",maxval(rhogq(:,:,I_QG))," min:",minval(rhogq(:,:,I_QG))
      write(IO_FID_LOG,*) "I_NC: max",maxval(rhogq(:,:,I_NC))," min:",minval(rhogq(:,:,I_NC))
      write(IO_FID_LOG,*) "I_NR: max",maxval(rhogq(:,:,I_NR))," min:",minval(rhogq(:,:,I_NR))
      write(IO_FID_LOG,*) "I_NI: max",maxval(rhogq(:,:,I_NI))," min:",minval(rhogq(:,:,I_NI))
      write(IO_FID_LOG,*) "I_NS: max",maxval(rhogq(:,:,I_NS))," min:",minval(rhogq(:,:,I_NS))
      write(IO_FID_LOG,*) "I_NG: max",maxval(rhogq(:,:,I_NG))," min:",minval(rhogq(:,:,I_NG))
      write(IO_FID_LOG,*) "tem : max",maxval(tem(:,:))       ," min:",minval(tem(:,:))
      write(IO_FID_LOG,*) "pre : max",maxval(pre(:,:))       ," min:",minval(pre(:,:)) ! 09/08/18 T.Mitsui
      write(IO_FID_LOG,*) "rho : max",maxval(rho(:,:))       ," min:",minval(rho(:,:)) ! 09/08/18 T.Mitsui
      write(IO_FID_LOG,*) "PLCccn: max",maxval(dt*PLCccn(:,:))," min:",minval(dt*PLCccn(:,:))
      write(IO_FID_LOG,*) "PNCccn: max",maxval(dt*PNCccn(:,:))," min:",minval(dt*PNCccn(:,:))
      write(IO_FID_LOG,*) "PLIccn: max",maxval(dt*PLIccn(:,:))," min:",minval(dt*PLIccn(:,:))
      write(IO_FID_LOG,*) "PNIccn: max",maxval(dt*PNIccn(:,:))," min:",minval(dt*PNIccn(:,:))

      return
    end subroutine debugreport_nucleation

    !---------------------------------------------------------------------------
    subroutine debugreport_phasechange
      implicit none

      write(IO_FID_LOG,*) "*** Phase Change/mod_mp_ndw6"

      write(IO_FID_LOG,*) "I_QV: max",maxval(rhogq(:,:,I_QV))," min:",minval(rhogq(:,:,I_QV))
      write(IO_FID_LOG,*) "I_QC: max",maxval(rhogq(:,:,I_QC))," min:",minval(rhogq(:,:,I_QC))
      write(IO_FID_LOG,*) "I_QR: max",maxval(rhogq(:,:,I_QR))," min:",minval(rhogq(:,:,I_QR))
      write(IO_FID_LOG,*) "I_QI: max",maxval(rhogq(:,:,I_QI))," min:",minval(rhogq(:,:,I_QI))
      write(IO_FID_LOG,*) "I_QS: max",maxval(rhogq(:,:,I_QS))," min:",minval(rhogq(:,:,I_QS))
      write(IO_FID_LOG,*) "I_QG: max",maxval(rhogq(:,:,I_QG))," min:",minval(rhogq(:,:,I_QG))

      write(IO_FID_LOG,*) "I_NC: max",maxval(rhogq(:,:,I_NC))," min:",minval(rhogq(:,:,I_NC))
      write(IO_FID_LOG,*) "I_NR: max",maxval(rhogq(:,:,I_NR))," min:",minval(rhogq(:,:,I_NR))
      write(IO_FID_LOG,*) "I_NI: max",maxval(rhogq(:,:,I_NI))," min:",minval(rhogq(:,:,I_NI))
      write(IO_FID_LOG,*) "I_NS: max",maxval(rhogq(:,:,I_NS))," min:",minval(rhogq(:,:,I_NS))
      write(IO_FID_LOG,*) "I_NG: max",maxval(rhogq(:,:,I_NG))," min:",minval(rhogq(:,:,I_NG))

      write(IO_FID_LOG,*) "tem : max",maxval(tem(:,:))       ," min:",minval(tem(:,:))
      write(IO_FID_LOG,*) "pre : max",maxval(pre(:,:))       ," min:",minval(pre(:,:)) ! 09/08/18 T.Mitsui
      write(IO_FID_LOG,*) "rho : max",maxval(rho(:,:))       ," min:",minval(rho(:,:)) ! 09/08/18 T.Mitsui

      write(IO_FID_LOG,*) "PLChom: max",maxval(wdt*1.D0*PLChom(:,:))," min:",minval(wdt*1.D0*PLChom(:,:))
      write(IO_FID_LOG,*) "PNChom: max",maxval(wdt*1.D0*PNChom(:,:))," min:",minval(wdt*1.D0*PNChom(:,:))
      write(IO_FID_LOG,*) "PLChet: max",maxval(wdt*1.D0*PLChet(:,:))," min:",minval(wdt*1.D0*PLChet(:,:))
      write(IO_FID_LOG,*) "PNChet: max",maxval(wdt*1.D0*PNChet(:,:))," min:",minval(wdt*1.D0*PNChet(:,:))
      write(IO_FID_LOG,*) "PLRhet: max",maxval(wdt*1.D0*PLRhet(:,:))," min:",minval(wdt*1.D0*PLRhet(:,:))
      write(IO_FID_LOG,*) "PNRhet: max",maxval(wdt*1.D0*PNRhet(:,:))," min:",minval(wdt*1.D0*PNRhet(:,:))

      write(IO_FID_LOG,*) "PLImlt: max",maxval(wdt*1.D0*PLImlt(:,:))," min:",minval(wdt*1.D0*PLImlt(:,:))
      write(IO_FID_LOG,*) "PNImlt: max",maxval(wdt*1.D0*PNImlt(:,:))," min:",minval(wdt*1.D0*PNImlt(:,:))
      write(IO_FID_LOG,*) "PLSmlt: max",maxval(wdt*1.D0*PLSmlt(:,:))," min:",minval(wdt*1.D0*PLSmlt(:,:))
      write(IO_FID_LOG,*) "PNSmlt: max",maxval(wdt*1.D0*PNSmlt(:,:))," min:",minval(wdt*1.D0*PNSmlt(:,:))
      write(IO_FID_LOG,*) "PLGmlt: max",maxval(wdt*1.D0*PLGmlt(:,:))," min:",minval(wdt*1.D0*PLGmlt(:,:))
      write(IO_FID_LOG,*) "PNGmlt: max",maxval(wdt*1.D0*PNGmlt(:,:))," min:",minval(wdt*1.D0*PNGmlt(:,:))

      write(IO_FID_LOG,*) "PLCdep: max",maxval(wdt*1.D0*PLCdep(:,:))," min:",minval(wdt*1.D0*PLCdep(:,:))
      write(IO_FID_LOG,*) "PLRdep: max",maxval(wdt*1.D0*PLRdep(:,:))," min:",minval(wdt*1.D0*PLRdep(:,:))
      write(IO_FID_LOG,*) "PLIdep: max",maxval(wdt*1.D0*PLIdep(:,:))," min:",minval(wdt*1.D0*PLIdep(:,:))
      write(IO_FID_LOG,*) "PLSdep: max",maxval(wdt*1.D0*PLSdep(:,:))," min:",minval(wdt*1.D0*PLSdep(:,:))
      write(IO_FID_LOG,*) "PLGdep: max",maxval(wdt*1.D0*PLGdep(:,:))," min:",minval(wdt*1.D0*PLGdep(:,:))

      write(IO_FID_LOG,*) "PNCdep: max",maxval(wdt*1.D0*PNCdep(:,:))," min:",minval(wdt*1.D0*PNCdep(:,:))
      write(IO_FID_LOG,*) "PNRdep: max",maxval(wdt*1.D0*PNRdep(:,:))," min:",minval(wdt*1.D0*PNRdep(:,:))
      write(IO_FID_LOG,*) "PNIdep: max",maxval(wdt*1.D0*PNIdep(:,:))," min:",minval(wdt*1.D0*PNIdep(:,:))
      write(IO_FID_LOG,*) "PNSdep: max",maxval(wdt*1.D0*PNSdep(:,:))," min:",minval(wdt*1.D0*PNSdep(:,:))
      write(IO_FID_LOG,*) "PNGdep: max",maxval(wdt*1.D0*PNGdep(:,:))," min:",minval(wdt*1.D0*PNGdep(:,:))

      return
    end subroutine debugreport_phasechange

    !---------------------------------------------------------------------------
    subroutine debugreport_collection
      implicit none

      write(IO_FID_LOG,*) "*** Collection process /mod_mp_ndw6"

      write(IO_FID_LOG,*) "I_QV: max",maxval(rhogq(:,:,I_QV))," min:",minval(rhogq(:,:,I_QV))
      write(IO_FID_LOG,*) "I_QC: max",maxval(rhogq(:,:,I_QC))," min:",minval(rhogq(:,:,I_QC))
      write(IO_FID_LOG,*) "I_QR: max",maxval(rhogq(:,:,I_QR))," min:",minval(rhogq(:,:,I_QR))
      write(IO_FID_LOG,*) "I_QI: max",maxval(rhogq(:,:,I_QI))," min:",minval(rhogq(:,:,I_QI))
      write(IO_FID_LOG,*) "I_QS: max",maxval(rhogq(:,:,I_QS))," min:",minval(rhogq(:,:,I_QS))
      write(IO_FID_LOG,*) "I_QG: max",maxval(rhogq(:,:,I_QG))," min:",minval(rhogq(:,:,I_QG))
      write(IO_FID_LOG,*) "I_NC: max",maxval(rhogq(:,:,I_NC))," min:",minval(rhogq(:,:,I_NC))
      write(IO_FID_LOG,*) "I_NR: max",maxval(rhogq(:,:,I_NR))," min:",minval(rhogq(:,:,I_NR))
      write(IO_FID_LOG,*) "I_NI: max",maxval(rhogq(:,:,I_NI))," min:",minval(rhogq(:,:,I_NI))
      write(IO_FID_LOG,*) "I_NS: max",maxval(rhogq(:,:,I_NS))," min:",minval(rhogq(:,:,I_NS))
      write(IO_FID_LOG,*) "I_NG: max",maxval(rhogq(:,:,I_NG))," min:",minval(rhogq(:,:,I_NG))

      write(IO_FID_LOG,*) "tem : max",maxval(tem(:,:))       ," min:",minval(tem(:,:))
      write(IO_FID_LOG,*) "pre : max",maxval(pre(:,:))       ," min:",minval(pre(:,:)) ! 09/08/18 T.Mitsui
      write(IO_FID_LOG,*) "rho : max",maxval(rho(:,:))       ," min:",minval(rho(:,:)) ! 09/08/18 T.Mitsui

      write(IO_FID_LOG,*) "PLCaut: max",maxval(wdt*1.D0*PLCaut(:,:))," min:",minval(wdt*1.D0*PLCaut(:,:))
      write(IO_FID_LOG,*) "PLCacc: max",maxval(wdt*1.D0*PLCacc(:,:))," min:",minval(wdt*1.D0*PLCacc(:,:))
      write(IO_FID_LOG,*) "PNCaut: max",maxval(wdt*1.D0*PNCaut(:,:))," min:",minval(wdt*1.D0*PNCaut(:,:))
      write(IO_FID_LOG,*) "PNCacc: max",maxval(wdt*1.D0*PNCacc(:,:))," min:",minval(wdt*1.D0*PNCacc(:,:))
      write(IO_FID_LOG,*) "PNRaut: max",maxval(wdt*1.D0*PNRaut(:,:))," min:",minval(wdt*1.D0*PNRaut(:,:))
      write(IO_FID_LOG,*) "PNRslc: max",maxval(wdt*1.D0*PNRslc(:,:))," min:",minval(wdt*1.D0*PNRslc(:,:))
      write(IO_FID_LOG,*) "PNRbrk: max",maxval(wdt*1.D0*PNRbrk(:,:))," min:",minval(wdt*1.D0*PNRbrk(:,:))
      ! cloud mass
      write(IO_FID_LOG,*) "PLGacLC2LG: max",maxval(wdt*1.D0*PLGacLC2LG(:,:))," min:",minval(wdt*1.D0*PLGacLC2LG(:,:))
      write(IO_FID_LOG,*) "PLSacLC2LS: max",maxval(wdt*1.D0*PLSacLC2LS(:,:))," min:",minval(wdt*1.D0*PLSacLC2LS(:,:))
      write(IO_FID_LOG,*) "PLIacLC2LI: max",maxval(wdt*1.D0*PLIacLC2LI(:,:))," min:",minval(wdt*1.D0*PLIacLC2LI(:,:))
      ! cloud num
      write(IO_FID_LOG,*) "PNGacNC2NG: max",maxval(wdt*1.D0*PNGacNC2NG(:,:))," min:",minval(wdt*1.D0*PNGacNC2NG(:,:))
      write(IO_FID_LOG,*) "PNIacNC2NI: max",maxval(wdt*1.D0*PNIacNC2NI(:,:))," min:",minval(wdt*1.D0*PNIacNC2NI(:,:))
      ! rain mass
      write(IO_FID_LOG,*) "PLRacLG2LG_R: max",maxval(wdt*1.D0*PLRacLG2LG(:,:))," min:",minval(wdt*1.D0*PLRacLG2LG(:,:))
      write(IO_FID_LOG,*) "PLRacLS2LG_S: max",maxval(wdt*1.D0*PLRacLS2LG_S(:,:))," min:",minval(wdt*1.D0*PLRacLS2LG_S(:,:))
      write(IO_FID_LOG,*) "PLRacLI2LG_R: max",maxval(wdt*1.D0*PLRacLI2LG_R(:,:))," min:",minval(wdt*1.D0*PLRacLI2LG_R(:,:))
      ! rain num
      write(IO_FID_LOG,*) "PNRacNI2NG_R: max",maxval(wdt*1.D0*PNRacNI2NG_R(:,:))," min:",minval(wdt*1.D0*PNRacNI2NG_R(:,:))
      write(IO_FID_LOG,*) "PNRacNS2NG_R: max",maxval(wdt*1.D0*PNRacNS2NG_R(:,:))," min:",minval(wdt*1.D0*PNRacNS2NG_R(:,:))
      write(IO_FID_LOG,*) "PNRacNG2NG  : max",maxval(wdt*1.D0*PNRacNG2NG  (:,:))," min:",minval(wdt*1.D0*PNRacNG2NG  (:,:))
      ! ice mass
      write(IO_FID_LOG,*) "PLRacLI2LG_I: max",maxval(wdt*1.D0*PLRacLI2LG_I(:,:))," min:",minval(wdt*1.D0*PLRacLI2LG_I(:,:))
      write(IO_FID_LOG,*) "PLIacLI2LS:   max",maxval(wdt*1.D0*PLIacLI2LS(:,:))," min:",minval(wdt*1.D0*PLIacLI2LS(:,:))
      write(IO_FID_LOG,*) "PLIacLS2LS_I: max",maxval(wdt*1.D0*PLIacLS2LS(:,:))," min:",minval(wdt*1.D0*PLIacLS2LS(:,:))
      ! ice num
      write(IO_FID_LOG,*) "PNRacNI2NG_I: max",maxval(wdt*1.D0*PNRacNI2NG_I(:,:))," min:",minval(wdt*1.D0*PNRacNI2NG_I(:,:))
      write(IO_FID_LOG,*) "PNIacNI2NS:   max",maxval(wdt*1.D0*PNIacNI2NS(:,:))," min:",minval(wdt*1.D0*PNIacNI2NS(:,:))
      write(IO_FID_LOG,*) "PNIacNS2NS_I: max",maxval(wdt*1.D0*PNIacNS2NS(:,:))," min:",minval(wdt*1.D0*PNIacNS2NS(:,:))
      ! snow mass
      write(IO_FID_LOG,*) "PNIacNS2NS_S: max",maxval(wdt*1.D0*PNIacNS2NS(:,:))," min:",minval(wdt*1.D0*PNIacNS2NS(:,:))
      write(IO_FID_LOG,*) "PLRacLS2LG_S: max",maxval(wdt*1.D0*PLRacLS2LG_S(:,:))," min:",minval(wdt*1.D0*PLRacLS2LG_S(:,:))
      write(IO_FID_LOG,*) "PLGacLS2LG:   max",maxval(wdt*1.D0*PLGacLS2LG(:,:)),  " min:",minval(wdt*1.D0*PLGacLS2LG(:,:))
      ! snow num
      write(IO_FID_LOG,*) "PNRacNS2NG_S: max",maxval(wdt*1.D0*PNRacNS2NG_S(:,:))," min:",minval(wdt*1.D0*PNRacNS2NG_S(:,:))
      write(IO_FID_LOG,*) "PNGacNS2NG:   max",maxval(wdt*1.D0*PNGacNS2NG(:,:))," min:",minval(wdt*1.D0*PNGacNS2NG(:,:))
      ! grauepl num 
      write(IO_FID_LOG,*) "PNRacNG2NG_G: max",maxval(wdt*1.D0*PNRacNG2NG(:,:))," min:",minval(wdt*1.D0*PNRacNG2NG(:,:))

      write(IO_FID_LOG,*) "PLIcon: max",maxval(wdt*1.D0*PLIcon(:,:))," min:",minval(wdt*1.D0*PLIcon(:,:))
      write(IO_FID_LOG,*) "PLScon: max",maxval(wdt*1.D0*PLScon(:,:))," min:",minval(wdt*1.D0*PLScon(:,:))
      write(IO_FID_LOG,*) "PNIcon: max",maxval(wdt*1.D0*PNIcon(:,:))," min:",minval(wdt*1.D0*PNIcon(:,:))
      write(IO_FID_LOG,*) "PNScon: max",maxval(wdt*1.D0*PNScon(:,:))," min:",minval(wdt*1.D0*PNScon(:,:))
       
      write(IO_FID_LOG,*) "PLIacm: max",maxval(wdt*1.D0*PLIacm(:,:))," min:",minval(wdt*1.D0*PLIacm(:,:))
      write(IO_FID_LOG,*) "PLIarm: max",maxval(wdt*1.D0*PLIarm(:,:))," min:",minval(wdt*1.D0*PLIarm(:,:))
      write(IO_FID_LOG,*) "PLSacm: max",maxval(wdt*1.D0*PLSacm(:,:))," min:",minval(wdt*1.D0*PLSacm(:,:))
      write(IO_FID_LOG,*) "PLSarm: max",maxval(wdt*1.D0*PLSarm(:,:))," min:",minval(wdt*1.D0*PLSarm(:,:))
      write(IO_FID_LOG,*) "PLGacm: max",maxval(wdt*1.D0*PLGacm(:,:))," min:",minval(wdt*1.D0*PLGacm(:,:))
      write(IO_FID_LOG,*) "PLGarm: max",maxval(wdt*1.D0*PLGarm(:,:))," min:",minval(wdt*1.D0*PLGarm(:,:))
       
      write(IO_FID_LOG,*) "PNIacm: max",maxval(wdt*1.D0*PNIacm(:,:))," min:",minval(wdt*1.D0*PNIacm(:,:))    
      write(IO_FID_LOG,*) "PNIarm: max",maxval(wdt*1.D0*PNIarm(:,:))," min:",minval(wdt*1.D0*PNIarm(:,:))    
      write(IO_FID_LOG,*) "PNSacm: max",maxval(wdt*1.D0*PNSacm(:,:))," min:",minval(wdt*1.D0*PNSacm(:,:))    
      write(IO_FID_LOG,*) "PNSarm: max",maxval(wdt*1.D0*PNSarm(:,:))," min:",minval(wdt*1.D0*PNSarm(:,:))    
      write(IO_FID_LOG,*) "PNGacm: max",maxval(wdt*1.D0*PNGacm(:,:))," min:",minval(wdt*1.D0*PNGacm(:,:))    
      write(IO_FID_LOG,*) "PNGarm: max",maxval(wdt*1.D0*PNGarm(:,:))," min:",minval(wdt*1.D0*PNGarm(:,:))    

      return
    end subroutine debugreport_collection

    !---------------------------------------------------------------------------
    subroutine debugreport_sedimentation
      implicit none
      !-------------------------------------------------------------------------

      write(IO_FID_LOG,*) "*** sedimentation /mod_mp_ndw6"

      write(IO_FID_LOG,*) "[Vt_QC] max:", maxval(velw(kmin:kmax,IS:IE,JS:JE,I_QC)), &
                                 " min:", minval(velw(kmin:kmax,IS:IE,JS:JE,I_QC))
      write(IO_FID_LOG,*) "[Vt_QR] max:", maxval(velw(kmin:kmax,IS:IE,JS:JE,I_QR)), &
                                 " min:", minval(velw(kmin:kmax,IS:IE,JS:JE,I_QR))
      write(IO_FID_LOG,*) "[Vt_QI] max:", maxval(velw(kmin:kmax,IS:IE,JS:JE,I_QI)), &
                                 " min:", minval(velw(kmin:kmax,IS:IE,JS:JE,I_QI))
      write(IO_FID_LOG,*) "[Vt_QS] max:", maxval(velw(kmin:kmax,IS:IE,JS:JE,I_QS)), &
                                 " min:", minval(velw(kmin:kmax,IS:IE,JS:JE,I_QS))
      write(IO_FID_LOG,*) "[Vt_QG] max:", maxval(velw(kmin:kmax,IS:IE,JS:JE,I_QG)), &
                                 " min:", minval(velw(kmin:kmax,IS:IE,JS:JE,I_QG))
      write(IO_FID_LOG,*) "[Vt_NC] max:", maxval(velw(kmin:kmax,IS:IE,JS:JE,I_NC)), &
                                 " min:", minval(velw(kmin:kmax,IS:IE,JS:JE,I_NC))
      write(IO_FID_LOG,*) "[Vt_NR] max:", maxval(velw(kmin:kmax,IS:IE,JS:JE,I_NR)), &
                                 " min:", minval(velw(kmin:kmax,IS:IE,JS:JE,I_NR))
      write(IO_FID_LOG,*) "[Vt_NI] max:", maxval(velw(kmin:kmax,IS:IE,JS:JE,I_NI)), &
                                 " min:", minval(velw(kmin:kmax,IS:IE,JS:JE,I_NI))
      write(IO_FID_LOG,*) "[Vt_NS] max:", maxval(velw(kmin:kmax,IS:IE,JS:JE,I_NS)), &
                                 " min:", minval(velw(kmin:kmax,IS:IE,JS:JE,I_NS))
      write(IO_FID_LOG,*) "[Vt_NG] max:", maxval(velw(kmin:kmax,IS:IE,JS:JE,I_NG)), &
                                 " min:", minval(velw(kmin:kmax,IS:IE,JS:JE,I_NG))

      write(IO_FID_LOG,*) "[Precip rain] max:", maxval(flux_rain(kmin-1,IS:IE,JS:JE)), &
                                       " min:", minval(flux_rain(kmin-1,IS:IE,JS:JE))
      write(IO_FID_LOG,*) "[Precip snow] max:", maxval(flux_rain(kmin-1,IS:IE,JS:JE)), &
                                       " min:", minval(flux_rain(kmin-1,IS:IE,JS:JE))

      return
    end subroutine debugreport_sedimentation

  end subroutine mp_ndw6

  !-----------------------------------------------------------------------------
  subroutine debug_tem( &
      point,    &
      tem,      &
      rho,      &
      pre,      &
      qv        )
    use mod_process, only: &
       PRC_myrank
    use mod_grid, only: &
       KA  => GRID_KA,  &
       KS  => GRID_KS,  &
       KE  => GRID_KE,  &
       IJA => GRID_IJA, &
       IJS => GRID_IJS, &
       IJE => GRID_IJE
    implicit none

    integer, intent(in) :: point
    real(8), intent(in) :: tem(IJA,KA)
    real(8), intent(in) :: rho(IJA,KA)
    real(8), intent(in) :: pre(IJA,KA)
    real(8), intent(in) :: qv (IJA,KA)

    integer :: k ,ij
    !---------------------------------------------------------------------------

    do k  = 1,  KA
    do ij = IJS, IJE
       if (      tem(ij,k) < tem_min &
            .OR. rho(ij,k) < rho_min &
            .OR. pre(ij,k) < 1.d0    ) then

          write(IO_FID_LOG,'(A,I3,A,4(F16.5),3(I6))') &
          "*** point: ", point, &
          " low tem,rho,pre:", tem(ij,k), rho(ij,k), pre(ij,k), qv(ij,k), ij, k, PRC_myrank
       endif
    enddo
    enddo

    return
  end subroutine debug_tem







  !-----------------------------------------------------------------------------
  subroutine mp_ndw6_diag_volume( &
       ijdim, kdim, kmin, kmax,&
       rho, &
       nc, nr, ni, ns, ng, &
       qc, qr, qi, qs, qg, &
       volc, volr, voli, vols, volg )
    use mod_atmos_cnst, only: &
         I_QC, I_QR, I_QI, I_QS, I_QG, &
         CNST_PI, &
         CNST_DWATR ! water density=1000[kg/m3]
    implicit none
    !
    integer, intent(in) :: ijdim
    integer, intent(in) :: kdim
    integer, intent(in) :: kmin
    integer, intent(in) :: kmax
    ! air density[kg/m3]
    real(8), intent(in) :: rho(ijdim,kdim)
    ! number concentration[/m3]
    real(8), intent(in) :: nc(ijdim,kdim)
    real(8), intent(in) :: nr(ijdim,kdim)
    real(8), intent(in) :: ni(ijdim,kdim)
    real(8), intent(in) :: ns(ijdim,kdim)
    real(8), intent(in) :: ng(ijdim,kdim)
    ! mixing ratio[kg/kg]
    real(8), intent(in) :: qc(ijdim,kdim)
    real(8), intent(in) :: qr(ijdim,kdim)
    real(8), intent(in) :: qi(ijdim,kdim)
    real(8), intent(in) :: qs(ijdim,kdim)
    real(8), intent(in) :: qg(ijdim,kdim)
    ! volume concentration[m3/m3]
    real(8), intent(out) :: volc(ijdim,kdim)
    real(8), intent(out) :: volr(ijdim,kdim)
    real(8), intent(out) :: voli(ijdim,kdim)
    real(8), intent(out) :: vols(ijdim,kdim)
    real(8), intent(out) :: volg(ijdim,kdim)
    !
    volc(:,:) = rho(:,:)*qc(:,:)/CNST_DWATR
    volr(:,:) = rho(:,:)*qr(:,:)/CNST_DWATR
    voli(:,:) = rho(:,:)*qi(:,:)/rhoi
    vols(:,:) = rho(:,:)*qs(:,:)/rhoi
    volg(:,:) = rho(:,:)*qg(:,:)/rhoi
    !
    return
  end subroutine mp_ndw6_diag_volume
  ! this is called from external procedure
  subroutine mp_ndw6_effective_radius(&
       ijdim, kdim, kmin, kmax,&
       rho, tem, pre,      &
       nc, nr, ni, ns, ng, &
       qc, qr, qi, qs, qg, &
       rec, rer, rei, res, reg )
    !
    use mod_atmos_cnst, only: &
         I_QC, I_QR, I_QI, I_QS, I_QG
    !
    implicit none
    !
    integer, intent(in) :: ijdim
    integer, intent(in) :: kdim
    integer, intent(in) :: kmin
    integer, intent(in) :: kmax
    ! atmospheric condition
    real(8), intent(in) :: rho(ijdim,kdim) ! air density[kg/m3]
    real(8), intent(in) :: tem(ijdim,kdim) ! temperature[K]
    real(8), intent(in) :: pre(ijdim,kdim) ! pressure[Pa]
    ! number concentration[/m3]
    real(8), intent(in) :: nc(ijdim,kdim)
    real(8), intent(in) :: nr(ijdim,kdim)
    real(8), intent(in) :: ni(ijdim,kdim)
    real(8), intent(in) :: ns(ijdim,kdim)
    real(8), intent(in) :: ng(ijdim,kdim)
    ! mixing ratio[kg/kg]
    real(8), intent(in) :: qc(ijdim,kdim)
    real(8), intent(in) :: qr(ijdim,kdim)
    real(8), intent(in) :: qi(ijdim,kdim)
    real(8), intent(in) :: qs(ijdim,kdim)
    real(8), intent(in) :: qg(ijdim,kdim)
    ! effective radius [m]
    real(8), intent(out) :: rec(ijdim,kdim)
    real(8), intent(out) :: rer(ijdim,kdim)
    real(8), intent(out) :: rei(ijdim,kdim)
    real(8), intent(out) :: res(ijdim,kdim)
    real(8), intent(out) :: reg(ijdim,kdim)
    !
    ! mass concentration[kg/m3] and mean particle mass[kg]
    real(8) :: lc(ijdim,kdim), xc(ijdim,kdim)
    real(8) :: lr(ijdim,kdim), xr(ijdim,kdim)
    real(8) :: li(ijdim,kdim), xi(ijdim,kdim)
    real(8) :: ls(ijdim,kdim), xs(ijdim,kdim)
    real(8) :: lg(ijdim,kdim), xg(ijdim,kdim)
    ! mean diameter[m]
    real(8) :: dc_xa(ijdim,kdim)     !
    real(8) :: dr_xa(ijdim,kdim)     !
    real(8) :: di_xa(ijdim,kdim)     !
    real(8) :: ds_xa(ijdim,kdim)     !
    real(8) :: dg_xa(ijdim,kdim)     !
    !
    real(8) :: wxr ! work
    real(8) :: ddr ! difference from equilibrium diameter 
    !
    real(8), parameter :: eps=1.d-30
    !
    integer :: ij, k
    !
    ! Preparation
    do k=1, kdim
       do ij=1, ijdim
          ! Mass concentration [kg/m3]
          lc(ij,k)     = rho(ij,k)*qc(ij,k) 
          lr(ij,k)     = rho(ij,k)*qr(ij,k) 
          li(ij,k)     = rho(ij,k)*qi(ij,k) 
          ls(ij,k)     = rho(ij,k)*qs(ij,k) 
          lg(ij,k)     = rho(ij,k)*qg(ij,k) 
          ! mean particle mass[kg]
          xc(ij,k)     = min(xc_max, max(xc_min, lc(ij,k)/(nc(ij,k)+nc_min) ))
          xr(ij,k)     = min(xr_max, max(xr_min, lr(ij,k)/(nr(ij,k)+nr_min) ))
          xi(ij,k)     = min(xi_max, max(xi_min, li(ij,k)/(ni(ij,k)+ni_min) ))
          xs(ij,k)     = min(xs_max, max(xs_min, ls(ij,k)/(ns(ij,k)+ns_min) ))
          xg(ij,k)     = min(xg_max, max(xg_min, lg(ij,k)/(ng(ij,k)+ng_min) ))
          ! diamter of average mass  SB06(32)
          dc_xa(ij,k)  = a_m(I_QC)*xc(ij,k)**b_m(I_QC)
          dr_xa(ij,k)  = a_m(I_QR)*xr(ij,k)**b_m(I_QR)
          di_xa(ij,k)  = a_m(I_QI)*xi(ij,k)**b_m(I_QI)
          ds_xa(ij,k)  = a_m(I_QS)*xs(ij,k)**b_m(I_QS)
          dg_xa(ij,k)  = a_m(I_QG)*xg(ij,k)**b_m(I_QG)
          !
       end do
    end do
    !
    call calc_effective_radius(        &
         ijdim, kdim, kmin, kmax,      & ! in
         tem,                          & ! in
         nc, nr, ni, ns, ng,           & ! in
         xi, xs, xg,                   & ! in ! [Add] 09/09/03 T.Mitsui
         dc_xa, dr_xa, di_xa, ds_xa, dg_xa, &
         rec, rer, rei, res, reg       )
    !
    return
  end subroutine mp_ndw6_effective_radius
  !
  subroutine calc_effective_radius( &
       ijdim, kdim, kmin, kmax,     &
       tem,                         &
       NC, NR, NI, NS, NG,          &
       xi, xs, xg,                  & ! in ! [Add] 09/09/03 T.Mitsui
       dc_ave, dr_ave, di_ave, ds_ave, dg_ave, &
!!$    re_liq, re_sol,              &
       re_qc, re_qr, re_qi, re_qs, re_qg )
    !
    use mod_atmos_cnst, only: &
         CNST_UNDEF, &
         CNST_PI,    &
         CNST_TEM00, &
         WLABEL,            &
         I_QC, I_QR, I_QI, I_QS, I_QG
    implicit none
    ! index
    integer, intent(in) :: ijdim
    integer, intent(in) :: kdim
    integer, intent(in) :: kmin
    integer, intent(in) :: kmax
    ! temperature[K]
    real(8), intent(in) :: tem(ijdim,kdim)
    ! number concentration[/m3]
    real(8), intent(in) :: NC(ijdim,kdim)
    real(8), intent(in) :: NR(ijdim,kdim)
    real(8), intent(in) :: NI(ijdim,kdim)
    real(8), intent(in) :: NS(ijdim,kdim)
    real(8), intent(in) :: NG(ijdim,kdim)
    !
    real(8), intent(in) :: xi(ijdim,kdim)
    real(8), intent(in) :: xs(ijdim,kdim)
    real(8), intent(in) :: xg(ijdim,kdim)
    ! diameter of average mass[kg/m3]
    real(8), intent(in) :: dc_ave(ijdim,kdim)
    real(8), intent(in) :: dr_ave(ijdim,kdim)
    real(8), intent(in) :: di_ave(ijdim,kdim)
    real(8), intent(in) :: ds_ave(ijdim,kdim)
    real(8), intent(in) :: dg_ave(ijdim,kdim)
    ! effective radius[m] of each hydrometeors
!!$    real(8), intent(out):: re_liq(ijdim,kdim) ! liquid all
!!$    real(8), intent(out):: re_sol(ijdim,kdim) ! solid all
    !
    real(8), intent(out):: re_qc(ijdim,kdim)  ! cloud only
    real(8), intent(out):: re_qr(ijdim,kdim)  ! rain only
    real(8), intent(out):: re_qi(ijdim,kdim)  ! ice only
    real(8), intent(out):: re_qs(ijdim,kdim)  ! snow only
    real(8), intent(out):: re_qg(ijdim,kdim)  ! graupel only
    !
    ! radius of average mass
    real(8) :: rc,rc2,rc3
    real(8) :: rr,rr2,rr3
    real(8) :: ri,ri2,ri3, di2,di3
    real(8) :: rs,rs2,rs3
    real(8) :: rg,rg2,rg3
    ! 2nd. and 3rd. order moment of DSD
    real(8) :: rc2m, rc3m
    real(8) :: rr2m, rr3m
    real(8) :: ri2m(ijdim,kdim), ri3m(ijdim,kdim)
    real(8) :: rs2m(ijdim,kdim), rs3m(ijdim,kdim)
    real(8) :: rg2m(ijdim,kdim), rg3m(ijdim,kdim)
    ! 
    real(8) :: r2m_solid
    real(8) :: r3m_solid
    ! work variables
    logical :: flag_rel(ijdim,kdim) ! liquic
    logical :: flag_rei(ijdim,kdim) ! ice
    real(8) :: work1,work2
    !
    real(8) :: coef_Fuetal1998=0.D0
    !
    real(8) :: r_pi 
    ! r2m_min is minimum value(moment of 1 particle with 1 micron)
    real(8), parameter :: r2m_min=1.d-12
    logical, save :: flag_first=.true.
    integer :: ij,k
    !
    r_pi = 1.d0/CNST_PI
    !
    re_qc(:,:)  = CNST_UNDEF
    re_qr(:,:)  = CNST_UNDEF
    re_qi(:,:)  = CNST_UNDEF
    re_qs(:,:)  = CNST_UNDEF
    re_qg(:,:)  = CNST_UNDEF
    !
!!$    re_liq(:,:) = CNST_UNDEF
!!$    re_sol(:,:) = CNST_UNDEF
    !
    ri2m(:,:) = 0.D0
    ri3m(:,:) = 0.D0
    rs2m(:,:) = 0.D0
    rs3m(:,:) = 0.D0
    rg2m(:,:) = 0.D0
    rg3m(:,:) = 0.D0
    !
    flag_rel(:,:)= .false.
    flag_rei(:,:)= .false.
    !
    !
    ! cloud, rain
    do k=kmin,kmax
       do ij=1, ijdim
          rc   = (0.5*dc_ave(ij,k))
          rr   = (0.5*dr_ave(ij,k))
          !
          rc2  = rc*rc
          rr2  = rr*rr
          rc3  = rc*rc2
          rr3  = rr*rr2
          !
          rc2m = coef_r2(I_QC)*rc2*NC(ij,k)
          rc3m = coef_r3(I_QC)*rc3*NC(ij,k)
          rr2m = coef_r2(I_QR)*rr2*NR(ij,k)
          rr3m = coef_r3(I_QR)*rr3*NR(ij,k)
          ! cloud effective radius
          if( rc >= rmin_re )then
             re_qc(ij,k)  = coef_re(I_QC)*rc
          end if
          ! rain effective radius
          if( rr >= rmin_re )then
             re_qr(ij,k)  = coef_re(I_QR)*rr
          end if
!!$          ! cloud + rain effective radius
!!$          if( rc2m+rr2m > r2m_min )then
!!$             re_liq(ij,k) = (rc3m+rr3m)/(rc2m+rr2m)
!!$             flag_rel(ij,k)=.true.
!!$          end if
       end do
    end do
    !
    ! Ac = pi*r*r 
    do k=kmin,kmax
       do ij=1, ijdim
          ri2m(ij,k)  = CNST_PI*coef_rea2(I_QI)*NI(ij,k)*a_rea2(I_QI)*xi(ij,k)**b_rea2(I_QI)
          rs2m(ij,k)  = CNST_PI*coef_rea2(I_QS)*NS(ij,k)*a_rea2(I_QS)*xs(ij,k)**b_rea2(I_QS)
          rg2m(ij,k)  = CNST_PI*coef_rea2(I_QG)*NG(ij,k)*a_rea2(I_QG)*xg(ij,k)**b_rea2(I_QG)
       end do
    end do
    !
    ! Fu(1996), eq.(3.11) or Fu etal.(1998), (2.5)
    coef_Fuetal1998 = 3.d0/(4.d0*rhoi)
    ri3m(:,:) = coef_Fuetal1998 * NI(:,:) * xi(:,:)
    rs3m(:,:) = coef_Fuetal1998 * NS(:,:) * xs(:,:)
    rg3m(:,:) = coef_Fuetal1998 * NG(:,:) * xg(:,:)
    !
    do k=kmin, kmax
       do ij=1, ijdim
          ! ice particles
          if( ri2m(ij,k) > r2m_min )then
             re_qi(ij,k) = ri3m(ij,k)/ri2m(ij,k)
          else
             re_qi(ij,k) = CNST_UNDEF
          end if
          ! snow particles
          if( rs2m(ij,k) > r2m_min )then
             re_qs(ij,k) = rs3m(ij,k)/rs2m(ij,k)
          else
             re_qs(ij,k) = CNST_UNDEF
          end if
          ! graupel particles
          if( rg2m(ij,k) > r2m_min )then
             re_qg(ij,k) = rg3m(ij,k)/rg2m(ij,k)
          else
             re_qg(ij,k) = CNST_UNDEF
          end if
          ! total solid particles
          r2m_solid = ri2m(ij,k)+rs2m(ij,k)+rg2m(ij,k)
          r3m_solid = ri3m(ij,k)+rs3m(ij,k)+rg3m(ij,k)
!!$          if( r2m_solid > r2m_min )then
!!$             re_sol(ij,k) = r3m_solid/r2m_solid
!!$             flag_rei(ij,k) = .true.
!!$          else
!!$             re_sol(ij,k) = CNST_UNDEF
!!$             flag_rei(ij,k) = .false.
!!$          end if
          !
       end do
    end do
    !
    return
  end subroutine calc_effective_radius
  !
  subroutine nucleation( &
       ijdim, kdim,      &
       kmin, kmax,       &
       l_region,         & 
       z, w,             & 
       rho, tem, pre,    & 
       LV, NC, NI,       &
       PNCccn, PLCccn,   &
       PNIccn, PLIccn,   &
       cva, cpa,         & ! in 
       dTdt_rad,         & ! in 
       qke,              & ! in 
       flag_history_in,  & ! in 
       ct, dt            ) ! in
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_atmos_cnst, only : &
         CNST_UNDEF,     & ! 10/08/03 T.Mitsui
         CNST_RVAP,      &
         CNST_EGRAV ! 09/08/18 T.Mitsui
    use mod_satadjust, only : &
         moist_psat_water,    &
         moist_psat_ice,      &
         moist_qsat_water,    & 
         moist_dqsw_dtem_rho, & 
         moist_dqsw_dtem_dpre,& 
         moist_qsat_ice,      & 
         moist_dqsi_dtem_rho    
    !
    implicit none
    !
    integer, intent(in)  :: ijdim
    integer, intent(in)  :: kdim
    integer, intent(in)  :: kmin
    integer, intent(in)  :: kmax
    integer, intent(in)  :: l_region ! 09/04/14x T.Mitsui
    ! 
    real(8), intent(in)  :: z(kdim)      ! 
    real(8), intent(in)  :: rho(ijdim,kdim)    ! [Add] 09/08/18 T.Mitsui
    real(8), intent(in)  :: tem(ijdim,kdim)    ! [Add] 09/08/18 T.Mitsui
    real(8), intent(in)  :: pre(ijdim,kdim)    ! [Add] 09/08/18 T.Mitsui
    real(8), intent(in)  :: w(ijdim,kdim)      ! w of half point
    !
    real(8), intent(in)  :: LV(ijdim,kdim)     !
    real(8), intent(in)  :: NC(ijdim,kdim)     ! [Add] 09/04/14 T.Mitsui
    real(8), intent(in)  :: NI(ijdim,kdim)     !
    real(8), intent(out) :: PNCccn(ijdim,kdim) !
    real(8), intent(out) :: PLCccn(ijdim,kdim) !
    real(8), intent(out) :: PNIccn(ijdim,kdim) !
    real(8), intent(out) :: PLIccn(ijdim,kdim)
    !
    real(8), intent(in) ::  cva(ijdim,kdim)      ! in  09/08/18 [Add] T.Mitsui
    real(8), intent(in) ::  cpa(ijdim,kdim)      ! in  09/08/18 [Add] T.Mitsui
    real(8), intent(in)  :: dTdt_rad(ijdim,kdim) ! 09/08/18 T.Mitsui
    real(8), intent(in)  :: qke(ijdim,kdim)      ! 09/08/18 T.Mitsui
    real(8), intent(in)  :: ct ! current time[sec.],   09/04/14x T.Mitsui
    real(8), intent(in)  :: dt
    logical, intent(in)  :: flag_history_in      ! in 10/08/03 [Add] T.Mitsui
    !
    ! namelist variables
    !
    ! total aerosol number concentration [/m3]
    real(8), parameter :: c_ccn_ocean= 1.00d8
    real(8), parameter :: c_ccn_land = 1.26d9
    real(8), save      :: c_ccn      = 1.00d8
    ! aerosol activation factor
    real(8), parameter :: kappa_ocean= 0.462d0
    real(8), parameter :: kappa_land = 0.308d0
    real(8), save      :: kappa      = 0.462d0
    real(8), save      :: c_in       = 1.d0
    ! SB06 (36)
    real(8), save :: nm_M92 = 1.d3
    real(8), save :: am_M92 = -0.639d0
    real(8), save :: bm_M92 = 12.96d0
    !
    real(8), save :: in_max = 1000.d3 ! max num. of Ice-Nuclei [num/m3]
    real(8), save :: ssi_max= 0.6d0
    real(8), save :: ssw_max= 1.1d0  ! [%]
    !
    logical, save :: flag_first = .true.
    real(8), save :: qke_min = 0.03d0 ! sigma=0.1[m/s], 09/08/18 T.Mitsui
    real(8), save :: tem_ccn_low=233.15d0  ! = -40 degC  ! [Add] 10/08/03 T.Mitsui
    real(8), save :: tem_in_low =173.15d0  ! = -100 degC ! [Add] 10/08/03 T.Mitsui
    !
    namelist /nm_mp_ndw6_nucleation/ &
         in_max,                     & ! 
         c_ccn, kappa,               & ! cloud nucleation
         nm_M92, am_M92, bm_M92,     & ! ice nucleation
         xc_ccn, xi_ccn,             &
         tem_ccn_low,                & ! [Add] 10/08/03 T.Mitsui
         tem_in_low,                 & ! [Add] 10/08/03 T.Mitsui
         ssw_max, ssi_max
    !
    real(8) :: c_ccn_map(ijdim,1)   ! c_ccn horizontal distribution
    real(8) :: kappa_map(ijdim,1)   ! kappa horizontal distribution
    real(8) :: c_in_map(ijdim,1)    ! c_in  horizontal distribution ! [Add] 11/08/30 T.Mitsui
    real(8) :: esw(ijdim,kdim)      ! saturation vapor pressure, water
    real(8) :: esi(ijdim,kdim)      !                            ice
    real(8) :: ssw(ijdim,kdim)      ! super saturation (water) 
    real(8) :: ssi(ijdim,kdim)      ! super saturation (ice) 
    real(8) :: w_dsswdz(ijdim,kdim) ! w*(d_ssw/ d_z) super saturation(water) flux
    real(8) :: w_dssidz(ijdim,kdim) ! w*(d_ssi/ d_z), 09/04/14 T.Mitsui
    real(8) :: ssw_below(ijdim,kdim)! ssw(k-1)
    real(8) :: ssi_below(ijdim,kdim)! ssi(k-1), 09/04/14 T.Mitsui
    real(8) :: z_below(ijdim,kdim)  ! z(k-1)
    real(8) :: dz                   ! z(k)-z(k-1)
    real(8) :: pv                   ! vapor pressure
    real(8) :: n_in                 ! number of ice nuclei
    ! work variables for Twomey Equation.
    real(8) :: qsw(ijdim,kdim)
    real(8) :: r_qsw(ijdim,kdim)
    real(8) :: qsi(ijdim,kdim)
    real(8) :: dqsidtem_rho(ijdim,kdim)
    real(8) :: dssidt_rad(ijdim,kdim)
    real(8) :: dni_ratio(ijdim,kdim)
    real(8) :: wssi, wdssi
    real(8) :: in0
    !
    real(8) :: xi_nuc(ijdim,1)    ! xi use the value @ cloud base
    real(8) :: alpha_nuc(ijdim,1) ! alpha_nuc 
    real(8) :: eta_nuc(ijdim,1)   ! xi use the value @ cloud base
    real(8) :: Dw, Q1, Q2
    !
    real(8) :: sigma_w(ijdim,kdim)
    real(8) :: weff(ijdim,kdim)
    real(8) :: weff_max(ijdim,kdim)
    real(8) :: w_m(ijdim,kdim)
    !    
    real(8) :: coef_ccn(ijdim)
    real(8) :: slope_ccn(ijdim)
    real(8) :: nc_new(ijdim,kdim)
    real(8) :: nc_new_below(ijdim,kdim)
    real(8) :: dnc_new
    real(8) :: nc_new_max   ! Lohmann (2002),
    real(8) :: a_max
    real(8) :: b_max
    logical :: flag_nucleation(ijdim,kdim)
    !
    real(8) :: r_gravity
    real(8), parameter :: r_sqrt3=0.577350269d0 ! = sqrt(1.d0/3.d0)    
    real(8), parameter :: eps=1.d-30 
    !====> ! 09/08/18
    !
    real(8) :: dlcdt_max, dli_max ! defined by supersaturation
    real(8) :: dncdt_max, dni_max ! defined by supersaturation
    real(8) :: rdt
    real(8) :: work
    !
    integer :: ij, k, kk
    !
    if( flag_first )then
       rewind(IO_FID_CONF)
       read(IO_FID_CONF, nml=nm_mp_ndw6_nucleation, end=100)
100    write(IO_FID_LOG, nml=nm_mp_ndw6_nucleation)
       flag_first=.false.
    end if
    !
    c_ccn_map(:,1) = c_ccn
    kappa_map(:,1) = kappa
    c_in_map(:,1)  = c_in
    !    
    nc_uplim(:,1)  = c_ccn_map(:,1)*1.5d0
    !
    rdt            = 1.d0/dt
    r_gravity      = 1.d0/CNST_EGRAV
    PNCccn(:,:)    = 0.D0
    PLCccn(:,:)    = 0.D0
    PNIccn(:,:)    = 0.D0
    PLIccn(:,:)    = 0.D0
    ssw(:,:)       = 0.D0
    ssi(:,:)       = 0.D0
    ssw_below(:,:) = 0.D0
    ssi_below(:,:) = 0.D0
    w_dsswdz(:,:)  = 0.D0
    w_dssidz(:,:)  = 0.D0
    dssidt_rad(:,:)= 0.D0
    dni_ratio(:,:) = CNST_UNDEF
    z_below(:,:)   = 0.D0
    weff(:,:)      = 0.D0
    work           = r_sqrt3*sqrt(qke_min)
    sigma_w(:,:)   = work
    nc_new(:,:)      = 0.D0
    nc_new_below(:,:)= 0.D0
    !      
    call moist_psat_water    ( tem, esw )
    call moist_psat_ice      ( tem, esi )
    call moist_qsat_water    ( tem, pre, qsw )
    call moist_qsat_ice      ( tem, pre, qsi ) 
    call moist_dqsi_dtem_rho ( tem, rho, dqsidtem_rho ) 
    r_qsw(:,:) = 1.d0/qsw
    !
    do k=kmin, kmax
       do ij=1, ijdim
          pv        = LV(ij,k)*CNST_RVAP*tem(ij,k)
          ssw(ij,k) = (pv/esw(ij,k) - 1.0d0)*100.D0
          ssi(ij,k) = (pv/esi(ij,k) - 1.0d0)
          ssw_below(ij,k+1) = ssw(ij,k)
          ssi_below(ij,k+1) = ssi(ij,k)
          z_below(ij,k+1)   = z(k)
       end do
    end do
    ssw_below(:,kmin) = ssw(:,kmin)
    ssi_below(:,kmin) = ssi(:,kmin)
    z_below(:,kmin)   = z(kmin-1)
    !
    !
    ! dS/dz is evaluated by first order upstream difference
    !***  Solution for Twomey Equation ***
    do ij=1, ijdim
       coef_ccn(ij)  = 1.d6*0.88d0*(c_ccn_map(ij,1)*1.d-6)**(2.d0/(kappa_map(ij,1) + 2.d0)) * &
            (70.D0)**(kappa_map(ij,1)/(kappa_map(ij,1) + 2.d0))
       slope_ccn(ij) = 1.5d0*kappa_map(ij,1)/(kappa_map(ij,1) + 2.d0)
    end do
    !
    do k=kmin, kmax
       sigma_w(:,k) = r_sqrt3*sqrt(max(qke(:,k),qke_min))
    end do
    sigma_w(:,kmin-1) = sigma_w(:,kmin)
    sigma_w(:,kmax+1) = sigma_w(:,kmax)
    ! effective vertical velocity 
    do k=kmin, kmax-1
       do ij=1, ijdim
          weff(ij,k) = 0.5d0*(w(ij,k) + w(ij,k+1)) - cpa(ij,k)*r_gravity*dTdt_rad(ij,k)
       end do
    end do
    weff(:,kmin-1) = weff(:,kmin)
    weff(:,kmax)   = weff(:,kmax-1)
    ! Lohmann (2002),JAS, eq.(1) but changing unit [cm-3] => [m-3]
    a_max = 1.d6*0.1d0*(1.d-6)**1.27d0
    b_max = 1.27d0 
    ! diagnose cloud condensation nuclei
    do k=kmin, kmax
       do ij=1, ijdim
          ! effective vertical velocity (maximum vertical velocity in turbulent flow) 
          weff_max(ij,k) = weff(ij,k) + sigma_w(ij,k) 
          ! large scale upward motion region and saturated
          if( (weff(ij,k) > 1.d-8) .and. (ssw(ij,k) > 1.d-10)  .and. pre(ij,k) > 300.d2 )then 
             ! Lohmann (2002), eq.(1)
             nc_new_max   = coef_ccn(ij)*weff_max(ij,k)**slope_ccn(ij)
             nc_new(ij,k) = a_max*nc_new_max**b_max 
          end if
       end do
    end do
    !
    flag_nucleation(:,:)=.false.
    do k=kmin, kmax
       do ij=1, ijdim
          ! nc_new is bound by upper limit
          if( nc_new(ij,k) > nc_uplim(ij,1) )then ! no more CCN
             flag_nucleation(ij,k) = .false.
             nc_new_below(ij,k+1)  = 1.d30
          else if( nc_new(ij,k) > eps )then ! nucleation can occur
             flag_nucleation(ij,k) = .true.
             nc_new_below(ij,k+1)  = nc_new(ij,k)
          else ! nucleation cannot occur(unsaturated or negative w)
             flag_nucleation(ij,k) = .false.
             nc_new_below(ij,k+1)  = 0.D0
          end if
       end do
    end do
    nc_new_below(:,kmin) = 0.D0 
    ! search maximum value of nc_new
    do k=kmin, kmax
       do ij=1, ijdim
          if(  ( nc_new(ij,k) < nc_new_below(ij,k) ) .or. &
               ( nc_new_below(ij,k) > c_ccn_map(ij,1)*0.05d0 ) )then ! 5% of c_ccn
             flag_nucleation(ij,k) = .false.
          end if
       end do
    end do
    ! nucleation occurs at only cloud base.
    ! if CCN is more than below parcel, nucleation newly occurs
    do k=kmin, kmax
       do ij=1, ijdim
          ! effective vertical velocity 
          if(   flag_nucleation(ij,k)               .and. & ! large scale upward motion region and saturated
               ( tem(ij,k)    > tem_ccn_low       ) .and. & 
               ( nc_new(ij,k) > NC(ij,k) )                )then
             dlcdt_max    = (LV(ij,k) - esw(ij,k)/(CNST_RVAP*tem(ij,k)))*rdt
             dncdt_max    = dlcdt_max/xc_min
             dnc_new      = nc_new(ij,k)-NC(ij,k)
             PNCccn(ij,k) = min( dncdt_max, dnc_new*rdt )
             PLCccn(ij,k) = min( dlcdt_max, xc_min*PNCccn(ij,k) )
          end if
       end do
    end do
    !
    ! ice nucleation
    !
    ssi_max = 1.d0
    ! +++ NOTE +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Based on Phillips etal.(2006).
    ! However this approach doesn't diagnose Ni itself but diagnose tendency.
    ! Original approach adjust Ni instantaneously .
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    do k=kmin+1, kmax
       do ij=1, ijdim
          dz             = z(k) - z_below(ij,k)
          w_dssidz(ij,k) = w(ij,k)*(ssi(ij,k) - ssi_below(ij,k))/dz ! 09/04/14 [Add] T.Mitsui
          dssidt_rad(ij,k) = -LV(ij,k)/(rho(ij,k)*qsi(ij,k)*qsi(ij,k))*dqsidtem_rho(ij,k)*dTdt_rad(ij,k) 
          dli_max        = (LV(ij,k) - esi(ij,k)/(CNST_RVAP*tem(ij,k)))*rdt
          dni_max        = min( dli_max/xi_ccn, (in_max-NI(ij,k))*rdt )
          wdssi          = min( w_dssidz(ij,k)+dssidt_rad(ij,k), 0.01d0) 
          wssi           = min( ssi(ij,k), ssi_max)                      
          ! SB06(34),(35)
          if(  ( wdssi    > eps        ) .and. & !
               (tem(ij,k) < 273.15d0   ) .and. & !
               (NI(ij,k)  < in_max     ) .and. &
               (wssi      >= eps       ) )then   !
             PNIccn(ij,k) = min(dni_max, c_in_map(ij,1)*bm_M92*nm_M92*0.3d0*exp(0.3d0*bm_M92*(wssi-0.1d0))*wdssi)
             PLIccn(ij,k) = min(dli_max, PNIccn(ij,k)*xi_ccn )
             ! only for output
             dni_ratio(ij,k) = dssidt_rad(ij,k)/( w_dssidz(ij,k)+dssidt_rad(ij,k) )
          else
             PNIccn(ij,k) = 0.D0
             PLIccn(ij,k) = 0.D0
          end if
          
       end do
    end do
    !
    return
  end subroutine nucleation
  ! 
  subroutine ice_multiplication( &
       ijdim, kdim,              & ! in
       kmin, kmax,               & ! in
       PLGspl, PLSspl, PNIspl,   & ! out
       PLIacLC2LI,               & ! in
       PLSacLC2LS,               & ! in
       PLGacLC2LG,               & ! in
       tem, LC, nc, xi, xc       ) ! in
    use mod_atmos_cnst, only: &
         I_QC
    !
    ! ice multiplication by splintering
    ! we consider Hallet-Mossop process
    implicit none
    !
    integer, intent(in) :: ijdim
    integer, intent(in) :: kdim
    integer, intent(in) :: kmin
    integer, intent(in) :: kmax
    real(8), intent(out):: PLGspl(ijdim,kdim) 
    real(8), intent(out):: PLSspl(ijdim,kdim) ! [Add]
    real(8), intent(out):: PNIspl(ijdim,kdim)
    !
    real(8), intent(in) :: PLIacLC2LI(ijdim, kdim)
    real(8), intent(in) :: PLSacLC2LS(ijdim, kdim)
    real(8), intent(in) :: PLGacLC2LG(ijdim, kdim)
    real(8), intent(in) :: tem(ijdim,kdim)
    real(8), intent(in) :: LC(ijdim,kdim)
    real(8), intent(in) :: NC(ijdim,kdim)
    real(8), intent(in) :: xi(ijdim,kdim)
    real(8), intent(in) :: xc(ijdim,kdim)
    !
    logical, save :: flag_first      = .true.
    ! production of (350.d3 ice particles)/(cloud riming [g]) => 350*d6 [/kg]
    real(8), parameter :: pice = 350.0d6
    ! production of (1 ice particle)/(250 cloud particles riming)
    real(8), parameter :: pnc  = 250.D0
    real(8), parameter :: rc_cr= 12.d-6 ! critical size[micron]
    ! temperature factor
    real(8) :: fp
    ! work for incomplete gamma function
    real(8), save :: xc_cr    ! mass[kg] of cloud with r=critical size[micron]
    real(8), save :: alpha    ! slope parameter of gamma function
    real(8), save :: gm, lgm  ! gamma(alpha), log(gamma(alpha))
    real(8) :: igm            ! in complete gamma(x,alpha)
    real(8) :: x
    ! coefficient of expansion using in calculation of igm
    real(8) :: a0,a1,a2,a3,a4,a5
    real(8) :: a6,a7,a8,a9,a10
    real(8) :: an1,an2,b0,b1,b2,c0,c1,c2
    real(8) :: d0,d1,d2,e0,e1,e2,h0,h1,h2
    real(8), parameter :: eps=1.0d-30
    ! number of cloud droplets larger than 12 micron(radius).
    real(8) :: n12            
    ! 
    real(8) :: wn, wni, wns, wng
    integer :: ij, k
    !
    if( flag_first )then
       flag_first = .false.
       ! work for Incomplete Gamma function
       xc_cr = (2.d0*rc_cr/a_m(I_QC))**(1.d0/b_m(I_QC))         
       alpha = (nu(I_QC)+1.d0)/mu(I_QC)
       gm    = gammafunc(alpha)
       lgm   = log(gm)
    end if
    PLGspl(:,1:kmin)   =0.D0
    PLSspl(:,1:kmin)   =0.D0
    PNIspl(:,1:kmin)   =0.D0
    PLGspl(:,kmax:kdim)=0.D0
    PLSspl(:,kmax:kdim)=0.D0
    PNIspl(:,kmax:kdim)=0.D0
    !
    do k=kmin, kmax
       do ij=1, ijdim
          ! Here we assume particle temperature is same as environment temperature.
          ! If you want to treat in a better manner, 
          ! you can diagnose with eq.(64) in CT(86)
          if     (tem(ij,k) >  270.16d0)then
             fp   = 0.D0
          else if(tem(ij,k) >= 268.16d0)then
             fp   = (270.16d0-tem(ij,k))*0.5d0
          else if(tem(ij,k) >= 265.16d0)then
             fp   = (tem(ij,k)-265.16d0)*0.333333333d0
          else
             fp   = 0.D0
          end if
          ! Approximation of Incomplete Gamma function
          ! Here we calculate with algorithm by Numerical Recipes.
          ! This approach is based on eq.(78) in Cotton etal.(1986), 
          ! but more accurate and expanded for Generalized Gamma Distribution.
          x       = coef_lambda(I_QC)*(xc_cr/xc(ij,k))**mu(I_QC)
          !
          if(x<1.d-2*alpha)then        ! negligible
             igm  = 0.D0
          else if(x<alpha+1.d0)then    ! series expansion
             ! 10th-truncation is enough for cloud droplet.
             a0   = 1.d0/alpha         ! n=0
             a1   = a0*x/(alpha+1.d0)  ! n=1
             a2   = a1*x/(alpha+2.d0)  ! n=2
             a3   = a2*x/(alpha+3.d0)  ! n=3
             a4   = a3*x/(alpha+4.d0)  ! n=4
             a5   = a4*x/(alpha+5.d0)  ! n=5
             a6   = a5*x/(alpha+6.d0)  ! n=6
             a7   = a6*x/(alpha+7.d0)  ! n=7
             a8   = a7*x/(alpha+8.d0)  ! n=8
             a9   = a8*x/(alpha+9.d0)  ! n=9
             a10  = a9*x/(alpha+10.D0) ! n=10
             igm  = (a0+a1+a2+a3+a4+a5+a6+a7+a8+a9+a10)*exp( -x + alpha*log(x) - lgm )
          else if(x<alpha*1.d2) then ! continued fraction expansion 
             ! 2nd-truncation is enough for cloud droplet.
             ! setup
             b0   = x+1.d0-alpha
             c0   = 1.d0/eps
             d0   = 1.d0/b0
             h0   = d0
             ! n=1
             an1  = -(1.d0-alpha)
             b1   = b0 + 2.d0
             d1   = 1.d0/(an1*d0+b1)
             c1   = b1+an1/c0
             e1   = d1*c1
             h1   = h0*e1
             ! n=2
             an2  = -2.d0*(2.d0-alpha)
             b2   = b1 + 2.d0
             d2   = 1.d0/(an2*d1+b2)
             c2   = b2+an2/c1
             e2   = d2*c2
             h2   = h1*e2
             !
             igm  = 1.d0 - exp( -x + alpha*log(x) - lgm )*h2
          else                       ! negligible
             igm  = 1.d0
          end if
          ! n12 is number of cloud droplets larger than 12 micron.
          n12 = NC(ij,k)*(1.d0-igm)
          ! eq.(82) CT(86)
          wn           = (pice + n12/((LC(ij,k)+xc_min)*pnc) )*fp ! filtered by xc_min
          wni          = wn*(-PLIacLC2LI(ij,k)  ) ! riming production rate is all negative
          wns          = wn*(-PLSacLC2LS(ij,k)  )
          wng          = wn*(-PLGacLC2LG(ij,k)  )
          PNIspl(ij,k) = wni+wns+wng
          !
          PLSspl(ij,k) = - wns*xi(ij,k) ! snow    => ice 
          PLGspl(ij,k) = - wng*xi(ij,k) ! graupel => ice
          !
       end do
    end do
    !
    return
  end subroutine ice_multiplication
  !
  subroutine mixed_phase_collection(   &
       ijdim, kdim,                    & ! in
       kmin, kmax,                     & ! in
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
       rho, & ! [Add] 11/08/30
       flag_history_in                 ) ! in [Add] 11/08/30
    use mod_stdio, only: &
       IO_FID_CONF,  &
       IO_FID_LOG,   &
       IO_L
    use mod_atmos_cnst, only: &
         pi => CNST_PI, &
         rhow => CNST_DWATR, &
         CNST_RVAP,     & 
         CNST_TEM00,    &
         CNST_CL,  &
         CNST_LHF0, &
         I_QC, I_QR, I_QI, I_QS, I_QG, &
         LHF 
    use mod_satadjust, only : &
         moist_psat_water,    &
         moist_psat_ice

    implicit none
    !
    integer, intent(in) :: ijdim
    integer, intent(in) :: kdim
    integer, intent(in) :: kmin
    integer, intent(in) :: kmax
    !--- mixed-phase collection process
    !    PXXacYY2ZZ_x: XX collecting YY to form ZZ.
    !                  _x means the contributions of XX or YY.
    !                  And all we set all production term as a negative sign to avoid confusion.
    ! 
    !--- ice-cloud     => ice
    real(8), intent(out):: PLIacLC2LI(ijdim,kdim) ! mass
    real(8), intent(out):: PNIacNC2NI(ijdim,kdim) ! number
    !--- snow-cloud    => snow
    real(8), intent(out):: PLSacLC2LS(ijdim,kdim) ! reduction of cloud
    real(8), intent(out):: PNSacNC2NS(ijdim,kdim) ! 
    !--- graupel-cloud => graupel
    real(8), intent(out):: PLGacLC2LG(ijdim,kdim)
    real(8), intent(out):: PNGacNC2NG(ijdim,kdim)
    !--- rain-ice      => graupel
    real(8), intent(out):: PLRacLI2LG_R(ijdim,kdim) ! reduction of rain
    real(8), intent(out):: PNRacNI2NG_R(ijdim,kdim) ! 
    real(8), intent(out):: PLRacLI2LG_I(ijdim,kdim) ! reduction of ice
    real(8), intent(out):: PNRacNI2NG_I(ijdim,kdim) ! 
    !--- rain-snow     => graupel
    real(8), intent(out):: PLRacLS2LG_R(ijdim,kdim) ! reduction of rain
    real(8), intent(out):: PNRacNS2NG_R(ijdim,kdim) ! 
    real(8), intent(out):: PLRacLS2LG_S(ijdim,kdim) ! reduction of snow
    real(8), intent(out):: PNRacNS2NG_S(ijdim,kdim) ! 
    !--- rain-graupel  => graupel
    real(8), intent(out):: PLRacLG2LG(ijdim,kdim) ! reduction of rain
    real(8), intent(out):: PNRacNG2NG(ijdim,kdim) ! reduction of graupel
    !--- ice-ice     => snow
    real(8), intent(out):: PLIacLI2LS(ijdim,kdim)
    real(8), intent(out):: PNIacNI2NS(ijdim,kdim)
    !--- ice-snow     => snow
    real(8), intent(out):: PLIacLS2LS(ijdim,kdim)
    real(8), intent(out):: PNIacNS2NS(ijdim,kdim)
    !--- snow-snow     => snow
    real(8), intent(out):: PNSacNS2NS(ijdim,kdim)
    !--- graupel-graupel=> graupel
    real(8), intent(out):: PNGacNG2NG(ijdim,kdim)
    !--- graupel-snow     => graupel
    real(8), intent(out):: PLGacLS2LG(ijdim,kdim)
    real(8), intent(out):: PNGacNS2NG(ijdim,kdim)
    !--- partial conversion
    !--- ice-cloud => graupel
    real(8), intent(out):: PLIcon(ijdim,kdim)
    real(8), intent(out):: PNIcon(ijdim,kdim)
    !--- snow-cloud => graupel
    real(8), intent(out):: PLScon(ijdim,kdim)
    real(8), intent(out):: PNScon(ijdim,kdim)
    !--- enhanced melting
    !--- graupel-cloud melting => rain
    real(8), intent(out):: PLGacm(ijdim,kdim)
    real(8), intent(out):: PNGacm(ijdim,kdim)
    !--- graupel-rain melting => rain
    real(8), intent(out):: PLGarm(ijdim,kdim)
    real(8), intent(out):: PNGarm(ijdim,kdim)
    !--- snow-cloud melting => rain
    real(8), intent(out):: PLSacm(ijdim,kdim)
    real(8), intent(out):: PNSacm(ijdim,kdim)
    !--- snow-rain melting => rain
    real(8), intent(out):: PLSarm(ijdim,kdim)
    real(8), intent(out):: PNSarm(ijdim,kdim)
    !--- ice-cloud melting => cloud ?
    real(8), intent(out):: PLIacm(ijdim,kdim)
    real(8), intent(out):: PNIacm(ijdim,kdim)
    !--- ice-rain melting => rain 
    real(8), intent(out):: PLIarm(ijdim,kdim)
    real(8), intent(out):: PNIarm(ijdim,kdim)
    ! 
    real(8), intent(in) :: wtem(ijdim,kdim)
    !--- mass concentration[kg/m3]
    real(8), intent(in) :: LV(ijdim,kdim)  
    real(8), intent(in) :: LC(ijdim,kdim)  
    real(8), intent(in) :: LR(ijdim,kdim)  
    real(8), intent(in) :: LI(ijdim,kdim)  
    real(8), intent(in) :: LS(ijdim,kdim)  
    real(8), intent(in) :: LG(ijdim,kdim)  
    !--- number concentration[/m3]
    real(8), intent(in) :: NC(ijdim,kdim)  
    real(8), intent(in) :: NR(ijdim,kdim)  
    real(8), intent(in) :: NI(ijdim,kdim)  
    real(8), intent(in) :: NS(ijdim,kdim)  
    real(8), intent(in) :: NG(ijdim,kdim)  
    ! necessary ?
    real(8), intent(in) :: xc(ijdim,kdim) ! LC/NC
    real(8), intent(in) :: xr(ijdim,kdim) ! LR/NR
    real(8), intent(in) :: xi(ijdim,kdim) ! LI/NI
    real(8), intent(in) :: xs(ijdim,kdim) ! LS/NS
    real(8), intent(in) :: xg(ijdim,kdim) ! LG/NG
    !--- diameter of averaged mass( D(ave x) )
    real(8), intent(in) :: dc_xave(ijdim,kdim)
    real(8), intent(in) :: dr_xave(ijdim,kdim)
    real(8), intent(in) :: di_xave(ijdim,kdim)
    real(8), intent(in) :: ds_xave(ijdim,kdim)
    real(8), intent(in) :: dg_xave(ijdim,kdim)
    !--- terminal velocity of averaged mass( vt(ave x) )
    real(8), intent(in) :: vtc_xave(ijdim,kdim)
    real(8), intent(in) :: vtr_xave(ijdim,kdim)
    real(8), intent(in) :: vti_xave(ijdim,kdim)
    real(8), intent(in) :: vts_xave(ijdim,kdim)
    real(8), intent(in) :: vtg_xave(ijdim,kdim)
    ! [Add] 11/08/30 T.Mitsui, for autoconversion of ice
    real(8), intent(in) :: rho(ijdim,kdim)
    logical, intent(in) :: flag_history_in
    !
    ! namelist variables 
    !=== for collection
    !--- threshold of diameters to collide with others
    real(8), save :: dc0 =  15.0d-6       ! lower threshold of cloud
    real(8), save :: dc1 =  40.0d-6       ! upper threshold of cloud
    real(8), save :: di0 = 150.0d-6       ! lower threshold of cloud
    real(8), save :: ds0 = 150.0d-6       ! lower threshold of cloud
    real(8), save :: dg0 = 150.0d-6       ! lower threshold of cloud
    !--- standard deviation of terminal velocity[m/s]
    real(8), save :: sigma_c=0.0d0        ! cloud
    real(8), save :: sigma_r=0.0d0        ! rain
    real(8), save :: sigma_i=0.2d0        ! ice  
    real(8), save :: sigma_s=0.2d0        ! snow
    real(8), save :: sigma_g=0.0d0        ! graupel
    !--- max collection efficiency for cloud
    real(8), save :: E_im = 0.80d0        ! ice max
    real(8), save :: E_sm = 0.80d0        ! snow max
    real(8), save :: E_gm = 1.00d0        ! graupel max
    !--- collection efficiency between 2 species
    real(8), save :: E_ir=1.d0            ! ice     x rain 
    real(8), save :: E_sr=1.d0            ! snow    x rain
    real(8), save :: E_gr=1.d0            ! graupel x rain
    real(8), save :: E_ii=1.d0            ! ice     x ice
    real(8), save :: E_si=1.d0            ! snow    x ice
    real(8), save :: E_gi=1.d0            ! graupel x ice
    real(8), save :: E_ss=1.d0            ! snow    x snow
    real(8), save :: E_gs=1.d0            ! graupel x snow
    real(8), save :: E_gg=1.d0            ! graupel x graupel
    !=== for partial conversion
    !--- flag: 1=> partial conversion to graupel, 0=> no conversion
    integer, save :: i_iconv2g=1          ! ice  => graupel 
    integer, save :: i_sconv2g=1          ! snow => graupel
    !--- bulk density of graupel 
    real(8), save :: rho_g   = 900.D0     ! [kg/m3]
    !--- space filling coefficient [%]
    real(8), save :: cfill_i = 0.68d0     ! ice 
    real(8), save :: cfill_s = 0.01d0     ! snow
    !--- critical diameter for ice conversion
    real(8), save :: di_cri  = 500.d-6    ! [m]
    ! [Add] 10/08/03 T.Mitsui
    logical, save :: opt_stick_KS96=.false.
    logical, save :: opt_stick_CO86=.false. 
    real(8), parameter :: a_dec = 0.883d0
    real(8), parameter :: b_dec = 0.093d0
    real(8), parameter :: c_dec = 0.00348d0
    real(8), parameter :: d_dec = 4.5185d-5
    !
    logical, save :: flag_first = .true.
    namelist /nm_mp_ndw6_collection/ &
         dc0, dc1, di0, ds0, dg0,    &
         sigma_c, sigma_r, sigma_i, sigma_s, sigma_g, &
         opt_stick_KS96,   &
         opt_stick_CO86,   & 
         E_im, E_sm, E_gm, &
         E_ir, E_sr, E_gr, E_ii, E_si, E_gi, E_ss, E_gs, E_gg, &
         i_iconv2g, i_sconv2g, rho_g, cfill_i, cfill_s, di_cri
    !
    real(8) :: tem(ijdim,kdim) 
    !
    !--- collection efficency of each specie
    real(8) :: E_c, E_r, E_i, E_s, E_g    ! 
    real(8) :: E_ic, E_sc, E_gc           !
    !--- sticking efficiency 
    real(8) :: E_stick(ijdim,kdim)
    ! [Add] 10/08/03 T.Mitsui
    real(8) :: temc, temc2, temc3
    real(8) :: E_dec
    real(8) :: esi_rat
    real(8) :: esi(ijdim,kdim)
    !
    real(8) :: temc_p, temc_m             ! celcius tem.
    ! [Add] 11/08/30 T.Mitsui, estimation of autoconversion time
    real(8) :: ci_aut(ijdim,kdim)
    real(8) :: taui_aut(ijdim,kdim)
    real(8) :: tau_sce(ijdim,kdim)
    !--- DSD averaged diameter for each species
    real(8) :: ave_dc                     ! cloud 
    real(8) :: ave_dr                     ! rain 
    real(8) :: ave_di                     ! ice 
    real(8) :: ave_ds                     ! snow
    real(8) :: ave_dg                     ! graupel
    !--- coefficient of collection equations(L:mass, N:number)
    real(8) :: coef_acc_LCI,   coef_acc_NCI   ! cloud     - cloud ice
    real(8) :: coef_acc_LCS,   coef_acc_NCS   ! cloud     - snow
    !
    real(8) :: coef_acc_LCG,   coef_acc_NCG   ! cloud     - graupel
    real(8) :: coef_acc_LRI_I, coef_acc_NRI_I ! rain      - cloud ice
    real(8) :: coef_acc_LRI_R, coef_acc_NRI_R ! rain      - cloud ice
    real(8) :: coef_acc_LRS_S, coef_acc_NRS_S ! rain      - snow
    real(8) :: coef_acc_LRS_R, coef_acc_NRS_R ! rain      - snow
    real(8) :: coef_acc_LRG,   coef_acc_NRG   ! rain      - graupel
    real(8) :: coef_acc_LII,   coef_acc_NII   ! cloud ice - cloud ice
    real(8) :: coef_acc_LIS,   coef_acc_NIS   ! cloud ice - snow
    real(8) ::                 coef_acc_NSS   ! snow      - snow
    real(8) ::                 coef_acc_NGG   ! grauepl   - graupel
    real(8) :: coef_acc_LSG,   coef_acc_NSG   ! snow      - graupel 
    !--- (diameter) x (diameter)
    real(8) :: dcdc, dcdi, dcds, dcdg
    real(8) :: drdr, drdi, drds, drdg
    real(8) :: didi, dids, didg
    real(8) :: dsds, dsdg
    real(8) :: dgdg
    !--- (terminal velocity) x (terminal velocity)
    real(8) :: vcvc, vcvi, vcvs, vcvg
    real(8) :: vrvr, vrvi, vrvs, vrvg
    real(8) :: vivi, vivs, vivg
    real(8) :: vsvs, vsvg
    real(8) :: vgvg
    !
    real(8) :: wx_cri, wx_crs
    real(8) :: coef_emelt
    real(8) :: w1
    !
    integer :: ij, k
    !
    if( flag_first )then
       rewind( IO_FID_CONF )
       read( IO_FID_CONF, nml=nm_mp_ndw6_collection, end=100 )
100    write( IO_FID_LOG, nml=nm_mp_ndw6_collection )
       flag_first = .false.
    end if
    !
    PLIacLC2LI(:,1:kmin)=0.D0
    PNIacNC2NI(:,1:kmin)=0.D0
    PLSacLC2LS(:,1:kmin)=0.D0
    PNSacNC2NS(:,1:kmin)=0.D0
    PLGacLC2LG(:,1:kmin)=0.D0
    PNGacNC2NG(:,1:kmin)=0.D0
    PLRacLI2LG_I(:,1:kmin)=0.D0
    PNRacNI2NG_I(:,1:kmin)=0.D0
    PLRacLI2LG_R(:,1:kmin)=0.D0
    PNRacNI2NG_R(:,1:kmin)=0.D0
    PLRacLS2LG_S(:,1:kmin)=0.D0
    PNRacNS2NG_S(:,1:kmin)=0.D0
    PLRacLS2LG_R(:,1:kmin)=0.D0
    PNRacNS2NG_R(:,1:kmin)=0.D0
    PLRacLG2LG(:,1:kmin)=0.D0
    PNRacNG2NG(:,1:kmin)=0.D0
    PLIacLI2LS(:,1:kmin)=0.D0
    PNIacNI2NS(:,1:kmin)=0.D0
    PLIacLS2LS(:,1:kmin)=0.D0
    PNIacNS2NS(:,1:kmin)=0.D0
    PNSacNS2NS(:,1:kmin)=0.D0
    PLGacLS2LG(:,1:kmin)=0.D0
    PNGacNS2NG(:,1:kmin)=0.D0
    PLIcon(:,1:kmin)=0.D0
    PNIcon(:,1:kmin)=0.D0
    PLScon(:,1:kmin)=0.D0
    PNScon(:,1:kmin)=0.D0
    PLIacm(:,1:kmin)=0.D0
    PNIacm(:,1:kmin)=0.D0
    PLIarm(:,1:kmin)=0.D0
    PNIarm(:,1:kmin)=0.D0
    PLSacm(:,1:kmin)=0.D0
    PNSacm(:,1:kmin)=0.D0
    PLSarm(:,1:kmin)=0.D0
    PNSarm(:,1:kmin)=0.D0
    PLGacm(:,1:kmin)=0.D0
    PNGacm(:,1:kmin)=0.D0
    PLGarm(:,1:kmin)=0.D0
    PNGarm(:,1:kmin)=0.D0
    PLIacLC2LI(:,kmax:kdim)=0.D0
    PNIacNC2NI(:,kmax:kdim)=0.D0
    PLSacLC2LS(:,kmax:kdim)=0.D0
    PNSacNC2NS(:,kmax:kdim)=0.D0
    PLGacLC2LG(:,kmax:kdim)=0.D0
    PNGacNC2NG(:,kmax:kdim)=0.D0
    PLRacLI2LG_I(:,kmax:kdim)=0.D0
    PNRacNI2NG_I(:,kmax:kdim)=0.D0
    PLRacLI2LG_R(:,kmax:kdim)=0.D0
    PNRacNI2NG_R(:,kmax:kdim)=0.D0
    PLRacLS2LG_S(:,kmax:kdim)=0.D0
    PNRacNS2NG_S(:,kmax:kdim)=0.D0
    PLRacLS2LG_R(:,kmax:kdim)=0.D0
    PNRacNS2NG_R(:,kmax:kdim)=0.D0
    PLRacLG2LG(:,kmax:kdim)=0.D0
    PNRacNG2NG(:,kmax:kdim)=0.D0
    PLIacLI2LS(:,kmax:kdim)=0.D0
    PNIacNI2NS(:,kmax:kdim)=0.D0
    PLIacLS2LS(:,kmax:kdim)=0.D0
    PNIacNS2NS(:,kmax:kdim)=0.D0
    PNSacNS2NS(:,kmax:kdim)=0.D0
    PLGacLS2LG(:,kmax:kdim)=0.D0
    PNGacNS2NG(:,kmax:kdim)=0.D0
    PLIcon(:,kmax:kdim)=0.D0
    PNIcon(:,kmax:kdim)=0.D0
    PLScon(:,kmax:kdim)=0.D0
    PNScon(:,kmax:kdim)=0.D0
    PLIacm(:,kmax:kdim)=0.D0
    PNIacm(:,kmax:kdim)=0.D0
    PLIarm(:,kmax:kdim)=0.D0
    PNIarm(:,kmax:kdim)=0.D0
    PLSacm(:,kmax:kdim)=0.D0
    PNSacm(:,kmax:kdim)=0.D0
    PLSarm(:,kmax:kdim)=0.D0
    PNSarm(:,kmax:kdim)=0.D0
    PLGacm(:,kmax:kdim)=0.D0
    PNGacm(:,kmax:kdim)=0.D0
    PLGarm(:,kmax:kdim)=0.D0
    PNGarm(:,kmax:kdim)=0.D0
    !
    ci_aut(:,:)   = 0.D0
    taui_aut(:,:) = 1.d10
    tau_sce(:,:)  = 1.d10
    !
    ! [Add] 10/08/03 T.Mitsui
    E_stick(:,:)=0.D0
    tem(:,:) = max(wtem(:,:), tem_min ) ! 11/08/30 T.Mitsui
    call moist_psat_ice( tem, esi )
    if( opt_stick_KS96 )then
       do k=kmin, kmax
          do ij=1, ijdim 
             ! Khain and Sednev (1996), eq.(3.15)
             temc          = tem(ij,k) - CNST_TEM00
             temc2         = temc*temc
             temc3         = temc2*temc
             E_dec         = max(0.D0, a_dec + b_dec*temc + c_dec*temc2 + d_dec*temc3 )
             esi_rat       = LV(ij,k)*CNST_RVAP*tem(ij,k)/esi(ij,k)
             E_stick(ij,k) = min(1.d0, E_dec*esi_rat)
          end do
       end do
    else if( opt_stick_CO86 )then 
       do k=kmin, kmax
          do ij=1, ijdim
             ! [Add] 11/08/30 T.Mitsui, Cotton et al. (1986)
             temc          = min(tem(ij,k) - CNST_TEM00,0.D0)
             w1            = 0.035d0*temc-0.7d0
             E_stick(ij,k) = 10.D0**w1
          end do
       end do
    else   
       do k=kmin, kmax
          do ij=1, ijdim 
             ! Lin et al. (1983)
             temc_m        = min(tem(ij,k) - CNST_TEM00,0.D0) ! T < 273.15
             E_stick(ij,k) = exp(0.09d0*temc_m)
          end do
       end do
    end if
    !
    do k=kmin, kmax
       do ij=1, ijdim 
          ! 
          temc_m = min(tem(ij,k) - CNST_TEM00,0.D0) ! T < 273.15
          temc_p = max(tem(ij,k) - CNST_TEM00,0.D0) ! T > 273.15
          ! averaged diameter using SB06(82)
          ave_dc = coef_d(I_QC)*xc(ij,k)**b_m(I_QC)
          ave_di = coef_d(I_QI)*xi(ij,k)**b_m(I_QI)
          ave_ds = coef_d(I_QS)*xs(ij,k)**b_m(I_QS)
          ave_dg = coef_d(I_QG)*xg(ij,k)**b_m(I_QG)
          !------------------------------------------------------------------------
          ! coellection efficiency are given as follows
          E_c = max(0.0d0, min(1.d0, (ave_dc-dc0)/(dc1-dc0) ))
          if(ave_di>di0)then
             E_i = E_im 
          else
             E_i = 0.D0
          end if
          if(ave_ds>ds0)then
             E_s = E_sm 
          else
             E_s = 0.D0
          end if
          if(ave_dg>dg0)then
             E_g = E_gm 
          else
             E_g = 0.D0
          end if
          E_ic = E_i*E_c 
          E_sc = E_s*E_c
          E_gc = E_g*E_c
          !------------------------------------------------------------------------
          ! Collection:  a collects b ( assuming particle size a>b )
          dcdc = dc_xave(ij,k) * dc_xave(ij,k) 
          drdr = dr_xave(ij,k) * dr_xave(ij,k) 
          didi = di_xave(ij,k) * di_xave(ij,k) 
          dsds = ds_xave(ij,k) * ds_xave(ij,k) 
          dgdg = dg_xave(ij,k) * dg_xave(ij,k) 
          dcdi = dc_xave(ij,k) * di_xave(ij,k) 
          dcds = dc_xave(ij,k) * ds_xave(ij,k) 
          dcdg = dc_xave(ij,k) * dg_xave(ij,k) 
          drdi = dr_xave(ij,k) * di_xave(ij,k) 
          drds = dr_xave(ij,k) * ds_xave(ij,k) 
          drdg = dr_xave(ij,k) * dg_xave(ij,k) 
          dids = di_xave(ij,k) * ds_xave(ij,k) 
          didg = di_xave(ij,k) * dg_xave(ij,k) 
          dsdg = ds_xave(ij,k) * dg_xave(ij,k) 
          !
          vcvc = vtc_xave(ij,k)* vtc_xave(ij,k)
          vrvr = vtr_xave(ij,k)* vtr_xave(ij,k)
          vivi = vti_xave(ij,k)* vti_xave(ij,k)
          vsvs = vts_xave(ij,k)* vts_xave(ij,k)
          vgvg = vtg_xave(ij,k)* vtg_xave(ij,k)
          vcvi = vtc_xave(ij,k)* vti_xave(ij,k)
          vcvs = vtc_xave(ij,k)* vts_xave(ij,k)
          vcvg = vtc_xave(ij,k)* vtg_xave(ij,k)
          vrvi = vtr_xave(ij,k)* vti_xave(ij,k)
          vrvs = vtr_xave(ij,k)* vts_xave(ij,k)
          vrvg = vtr_xave(ij,k)* vtg_xave(ij,k)
          vivs = vti_xave(ij,k)* vts_xave(ij,k)
          vivg = vti_xave(ij,k)* vtg_xave(ij,k)
          vsvg = vts_xave(ij,k)* vtg_xave(ij,k)
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
               +  sigma_i + sigma_c )**0.5d0
          coef_acc_NCI = &
               (  delta_b0(I_QC)*dcdc + delta_ab0(I_QI,I_QC)*dcdi + delta_b0(I_QI)*didi ) &
               *( theta_b0(I_QC)*vcvc - theta_ab0(I_QI,I_QC)*vcvi + theta_b0(I_QI)*vivi   &
               +  sigma_i + sigma_c )**0.5d0
          PLIacLC2LI(ij,k)= -0.25d0*pi*E_ic*NI(ij,k)*LC(ij,k)*coef_acc_LCI
          PNIacNC2NI(ij,k)= -0.25d0*pi*E_ic*NI(ij,k)*NC(ij,k)*coef_acc_NCI          
          ! cloud-snow => snow
          ! reduction term of cloud
          coef_acc_LCS = &
               (  delta_b1(I_QC)*dcdc + delta_ab1(I_QS,I_QC)*dcds + delta_b0(I_QS)*dsds ) &
               *( theta_b1(I_QC)*vcvc - theta_ab1(I_QS,I_QC)*vcvs + theta_b0(I_QS)*vsvs   &
               +  sigma_s + sigma_c )**0.5d0
          coef_acc_NCS = &
               (  delta_b0(I_QC)*dcdc + delta_ab0(I_QS,I_QC)*dcds + delta_b0(I_QS)*dsds ) &
               *( theta_b0(I_QC)*vcvc - theta_ab0(I_QS,I_QC)*vcvs + theta_b0(I_QS)*vsvs   &
               +  sigma_s + sigma_c )**0.5d0
          PLSacLC2LS(ij,k)= -0.25d0*pi*E_sc*NS(ij,k)*LC(ij,k)*coef_acc_LCS
          PNSacNC2NS(ij,k)= -0.25d0*pi*E_sc*NS(ij,k)*NC(ij,k)*coef_acc_NCS
          ! cloud-graupel => graupel
          ! reduction term of cloud
          coef_acc_LCG = &
               (  delta_b1(I_QC)*dcdc + delta_ab1(I_QG,I_QC)*dcdg + delta_b0(I_QG)*dgdg ) &
               *( theta_b1(I_QC)*vcvc - theta_ab1(I_QG,I_QC)*vcvg + theta_b0(I_QG)*vgvg   &
               +  sigma_g + sigma_c )**0.5d0
          coef_acc_NCG = &
               (  delta_b0(I_QC)*dcdc + delta_ab0(I_QG,I_QC)*dcdg + delta_b0(I_QG)*dgdg ) &
               *( theta_b0(I_QC)*vcvc - theta_ab0(I_QG,I_QC)*vcvg + theta_b0(I_QG)*vgvg   &
               +  sigma_g + sigma_c )**0.5d0
          PLGacLC2LG(ij,k)= -0.25d0*pi*E_gc*NG(ij,k)*LC(ij,k)*coef_acc_LCG
          PNGacNC2NG(ij,k)= -0.25d0*pi*E_gc*NG(ij,k)*NC(ij,k)*coef_acc_NCG
          ! snow-graupel => graupel
          coef_acc_LSG = &
               (  delta_b1(I_QS)*dsds + delta_ab1(I_QG,I_QS)*dsdg + delta_b0(I_QG)*dgdg ) &
               *( theta_b1(I_QS)*vsvs - theta_ab1(I_QG,I_QS)*vsvg + theta_b0(I_QG)*vgvg   &
               +  sigma_g + sigma_s )**0.5d0
          coef_acc_NSG = &
               (  delta_b0(I_QS)*dsds + delta_ab0(I_QG,I_QS)*dsdg + delta_b0(I_QG)*dgdg ) &
               ! [fix] T.Mitsui 08/05/08
               *( theta_b0(I_QS)*vsvs - theta_ab0(I_QG,I_QS)*vsvg + theta_b0(I_QG)*vgvg   &
               +  sigma_g + sigma_s )**0.5d0
          PLGacLS2LG(ij,k)= -0.25d0*pi*E_stick(ij,k)*E_gs*NG(ij,k)*LS(ij,k)*coef_acc_LSG
          PNGacNS2NG(ij,k)= -0.25d0*pi*E_stick(ij,k)*E_gs*NG(ij,k)*NS(ij,k)*coef_acc_NSG
          !------------------------------------------------------------------------
          ! ice-snow => snow
          ! reduction term of ice
          coef_acc_LIS = &
               (  delta_b1(I_QI)*didi + delta_ab1(I_QS,I_QI)*dids + delta_b0(I_QS)*dsds ) &
               *( theta_b1(I_QI)*vivi - theta_ab1(I_QS,I_QI)*vivs + theta_b0(I_QS)*vsvs   &
               +  sigma_i + sigma_s )**0.5d0
          coef_acc_NIS = &
               (  delta_b0(I_QI)*didi + delta_ab0(I_QS,I_QI)*dids + delta_b0(I_QS)*dsds ) &
               *( theta_b0(I_QI)*vivi - theta_ab0(I_QS,I_QI)*vivs + theta_b0(I_QS)*vsvs   &
               +  sigma_i + sigma_s )**0.5d0
          PLIacLS2LS(ij,k)= -0.25d0*pi*E_stick(ij,k)*E_si*NS(ij,k)*LI(ij,k)*coef_acc_LIS
          PNIacNS2NS(ij,k)= -0.25d0*pi*E_stick(ij,k)*E_si*NS(ij,k)*NI(ij,k)*coef_acc_NIS
          ! 
          if ( tem(ij,k) <= CNST_TEM00 )then
             ! rain-graupel => graupel
             ! reduction term of rain
             coef_acc_LRG = &
                  (  delta_b1(I_QR)*drdr + delta_ab1(I_QG,I_QR)*drdg + delta_b0(I_QG)*dgdg ) &
                  *( theta_b1(I_QR)*vrvr - theta_ab1(I_QG,I_QR)*vrvg + theta_b0(I_QG)*vgvg   &
                  +  sigma_r + sigma_g )**0.5d0
             PLRacLG2LG(ij,k) = -0.25d0*pi*E_gr*NG(ij,k)*LR(ij,k)*coef_acc_LRG
          else
             ! rain-graupel => rain
             ! reduction term of graupel
             coef_acc_LRG = &
                  (  delta_b1(I_QG)*dgdg + delta_ab1(I_QR,I_QG)*drdg + delta_b0(I_QR)*drdr ) &
                  *( theta_b1(I_QG)*vgvg - theta_ab1(I_QR,I_QG)*vrvg + theta_b0(I_QR)*vrvr   &
                  +  sigma_r + sigma_g )**0.5d0
             PLRacLG2LG(ij,k) = -0.25d0*pi*E_gr*NR(ij,k)*LG(ij,k)*coef_acc_LRG
          end if
          coef_acc_NRG = &
               (  delta_b0(I_QR)*drdr + delta_ab0(I_QG,I_QR)*drdg + delta_b0(I_QG)*dgdg ) &
               *( theta_b0(I_QR)*vrvr - theta_ab0(I_QG,I_QR)*vrvg + theta_b0(I_QG)*vgvg   &
               +  sigma_r + sigma_g )**0.5d0
          PNRacNG2NG(ij,k) = -0.25d0*pi*E_gr*NG(ij,k)*NR(ij,k)*coef_acc_NRG
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
               +  sigma_r + sigma_i )**0.5d0
          coef_acc_NRI_I = &
               (  delta_b0(I_QI)*didi + delta_ab0(I_QR,I_QI)*drdi + delta_b0(I_QR)*drdr ) &
               *( theta_b0(I_QI)*vivi - theta_ab0(I_QR,I_QI)*vrvi + theta_b0(I_QR)*vrvr   &
               +  sigma_r + sigma_i )**0.5d0
          PLRacLI2LG_I(ij,k)= -0.25d0*pi*E_ir*NR(ij,k)*LI(ij,k)*coef_acc_LRI_I
          PNRacNI2NG_I(ij,k)= -0.25d0*pi*E_ir*NR(ij,k)*NI(ij,k)*coef_acc_NRI_I
          ! reduction term of rain
          coef_acc_LRI_R = &
               (  delta_b1(I_QR)*drdr + delta_ab1(I_QI,I_QR)*drdi + delta_b0(I_QI)*didi ) &
               *( theta_b1(I_QR)*vrvr - theta_ab1(I_QI,I_QR)*vrvi + theta_b0(I_QI)*vivi   &
               +  sigma_r + sigma_i )**0.5d0
          coef_acc_NRI_R = &
               (  delta_b0(I_QR)*drdr + delta_ab0(I_QI,I_QR)*drdi + delta_b0(I_QI)*didi ) &
               *( theta_b0(I_QR)*vrvr - theta_ab0(I_QI,I_QR)*vrvi + theta_b0(I_QI)*vivi   &
               +  sigma_r + sigma_i )**0.5d0
          PLRacLI2LG_R(ij,k)= -0.25d0*pi*E_ir*NI(ij,k)*LR(ij,k)*coef_acc_LRI_R
          PNRacNI2NG_R(ij,k)= -0.25d0*pi*E_ir*NI(ij,k)*NR(ij,k)*coef_acc_NRI_R
          ! rain-snow => graupel
          ! reduction term of snow
          coef_acc_LRS_S = &
               (  delta_b1(I_QS)*dsds + delta_ab1(I_QR,I_QS)*drds + delta_b0(I_QR)*drdr ) &
               *( theta_b1(I_QS)*vsvs - theta_ab1(I_QR,I_QS)*vrvs + theta_b0(I_QR)*vrvr   &
               +  sigma_r + sigma_s )**0.5d0
          coef_acc_NRS_S = &
               (  delta_b0(I_QS)*dsds + delta_ab0(I_QR,I_QS)*drds + delta_b0(I_QR)*drdr ) &
               *( theta_b0(I_QS)*vsvs - theta_ab0(I_QR,I_QS)*vrvs + theta_b0(I_QR)*vrvr   &
               +  sigma_r + sigma_s )**0.5d0
          PLRacLS2LG_S(ij,k)= -0.25d0*pi*E_sr*NR(ij,k)*LS(ij,k)*coef_acc_LRS_S
          PNRacNS2NG_S(ij,k)= -0.25d0*pi*E_sr*NR(ij,k)*NS(ij,k)*coef_acc_NRS_S
          ! reduction term of rain
          coef_acc_LRS_R = &
               (  delta_b1(I_QR)*drdr + delta_ab1(I_QS,I_QR)*drds + delta_b0(I_QS)*dsds ) &
               *( theta_b1(I_QR)*vrvr - theta_ab1(I_QS,I_QR)*vrvs + theta_b0(I_QS)*vsvs   &
               +  sigma_r + sigma_s )**0.5d0
          coef_acc_NRS_R = &
               (  delta_b0(I_QR)*drdr + delta_ab0(I_QS,I_QR)*drds + delta_b0(I_QS)*dsds ) &
               *( theta_b0(I_QR)*vrvr - theta_ab0(I_QS,I_QR)*vrvs + theta_b0(I_QS)*vsvs   &
               +  sigma_r + sigma_s )**0.5d0
          PLRacLS2LG_R(ij,k)= -0.25d0*pi*E_sr*NS(ij,k)*LR(ij,k)*coef_acc_LRS_R
          PNRacNS2NG_R(ij,k)= -0.25d0*pi*E_sr*NS(ij,k)*NR(ij,k)*coef_acc_NRS_R
          !------------------------------------------------------------------------
          !
          !+++ pattern 3: a + a => b  (i-i)
          !                           
          !------------------------------------------------------------------------
          ! ice-ice ( reduction is double, but includes double-count)
          coef_acc_LII = &
               (  delta_b0(I_QI)*didi + delta_ab1(I_QI,I_QI)*didi + delta_b1(I_QI)*didi ) &
               *( theta_b0(I_QI)*vivi - theta_ab1(I_QI,I_QI)*vivi + theta_b1(I_QI)*vivi   &
               +  sigma_i + sigma_i )**0.5d0
          coef_acc_NII = &
               (  delta_b0(I_QI)*didi + delta_ab0(I_QI,I_QI)*didi + delta_b0(I_QI)*didi ) &
               *( theta_b0(I_QI)*vivi - theta_ab0(I_QI,I_QI)*vivi + theta_b0(I_QI)*vivi   &
               +  sigma_i + sigma_i )**0.5d0
          PLIacLI2LS(ij,k)= -0.25d0*pi*E_stick(ij,k)*E_ii*NI(ij,k)*LI(ij,k)*coef_acc_LII
          PNIacNI2NS(ij,k)= -0.25d0*pi*E_stick(ij,k)*E_ii*NI(ij,k)*NI(ij,k)*coef_acc_NII
          !
          ci_aut(ij,k)   =  0.25d0*pi*E_ii*NI(ij,k)*coef_acc_LII
          taui_aut(ij,k) = 1.d0/max(E_stick(ij,k)*ci_aut(ij,k),1.d-10)
          tau_sce(ij,k)  = LI(ij,k)/max(LI(ij,k)+LS(ij,k),1.d-10)
          !------------------------------------------------------------------------
          !
          !+++ pattern 4: a + a => a  (s-s)
          !
          !------------------------------------------------------------------------
          ! snow-snow => snow
          coef_acc_NSS = &
               (  delta_b0(I_QS)*dsds + delta_ab0(I_QS,I_QS)*dsds + delta_b0(I_QS)*dsds ) &
               *( theta_b0(I_QS)*vsvs - theta_ab0(I_QS,I_QS)*vsvs + theta_b0(I_QS)*vsvs   &
               +  sigma_s + sigma_s )**0.5d0
          PNSacNS2NS(ij,k)= -0.125d0*pi*E_stick(ij,k)*E_ss*NS(ij,k)*NS(ij,k)*coef_acc_NSS
          !
          ! graupel-grauple => graupel
          coef_acc_NGG = &
               (  delta_b0(I_QG)*dgdg + delta_ab0(I_QG,I_QG)*dgdg + delta_b0(I_QG)*dgdg ) &
               *( theta_b0(I_QG)*vgvg - theta_ab0(I_QG,I_QG)*vgvg + theta_b0(I_QG)*vgvg   &
               +  sigma_g + sigma_g )**0.5d0
          PNGacNG2NG(ij,k)= -0.125d0*pi*E_stick(ij,k)*E_gg*NG(ij,k)*NG(ij,k)*coef_acc_NGG
          !
          !------------------------------------------------------------------------
          !--- Partial conversion          
          ! SB06(70),(71)
          ! i_iconv2g: option whether partial conversions work or not
          ! ice-cloud => graupel
          if( ave_di > di_cri )then
             wx_cri = cfill_i*rhow/rho_g*( pi/6.d0*rho_g*ave_di*ave_di*ave_di/xi(ij,k) - 1.d0 ) 
             PLIcon(ij,k) = i_iconv2g*  PLIacLC2LI(ij,k)/max(1.d0, wx_cri)
             PNIcon(ij,k) = i_iconv2g*  PLIcon(ij,k)/xi(ij,k)      
          else
             wx_cri       = 0.D0
             PLIcon(ij,k) = 0.D0
             PNIcon(ij,k) = 0.D0
          end if
          ! snow-cloud => graupel
          wx_crs = cfill_s*rhow/rho_g*( pi/6.d0*rho_g*ave_ds*ave_ds*ave_ds/xs(ij,k) - 1.d0 ) 
          PLScon(ij,k) = i_sconv2g*  (PLSacLC2LS(ij,k))/max(1.d0, wx_crs)
          PNScon(ij,k) = i_sconv2g*  PLScon(ij,k)/xs(ij,k)
          !------------------------------------------------------------------------
          !--- enhanced melting( due to collection-freezing of water droplets )
          !    originally from Rutledge and Hobbs(1984). eq.(A.21)
          ! if T > 273.15 then temc_p=T-273.15, else temc_p=0
          ! 08/05/08 [fix] T.Mitsui LHF00 => LHF0
          ! melting occurs around T=273K, so LHF0 is suitable both SIMPLE and EXACT,
          ! otherwise LHF can have sign both negative(EXACT) and positive(SIMPLE).
!!$       coef_emelt   = -CNST_CL/CNST_LHF00*temc_p
          coef_emelt   =  CNST_CL/CNST_LHF0*temc_p
          ! cloud-graupel
          PLGacm(ij,k) =  coef_emelt*PLGacLC2LG(ij,k) 
          PNGacm(ij,k) =  PLGacm(ij,k)/xg(ij,k)
          ! rain-graupel
          PLGarm(ij,k) =  coef_emelt*PLRacLG2LG(ij,k)
          PNGarm(ij,k) =  PLGarm(ij,k)/xg(ij,k)
          ! cloud-snow
          PLSacm(ij,k) =  coef_emelt*(PLSacLC2LS(ij,k))
          PNSacm(ij,k) =  PLSacm(ij,k)/xs(ij,k)
          ! rain-snow
          PLSarm(ij,k) =  coef_emelt*(PLRacLS2LG_R(ij,k)+PLRacLS2LG_S(ij,k))
          PNSarm(ij,k) =  PLSarm(ij,k)/xg(ij,k)
          ! cloud-ice
          PLIacm(ij,k) =  coef_emelt*PLIacLC2LI(ij,k)
          PNIacm(ij,k) =  PLIacm(ij,k)/xi(ij,k)
          ! rain-ice
          PLIarm(ij,k) =  coef_emelt*(PLRacLI2LG_R(ij,k)+PLRacLI2LG_I(ij,k))
          PNIarm(ij,k) =  PLIarm(ij,k)/xg(ij,k)
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
  end subroutine mixed_phase_collection
  ! Auto-conversion, Accretion, Self-collection, Break-up
  subroutine aut_acc_slc_brk(  &
       ijdim, kdim,            &
       kmin, kmax,             &
       PLCaut, PNCaut,         &
       PNRaut,                 &
       PLCacc, PNCacc,         &
       PNRslc,                 &
       PNRbrk,                 &
       LC, LR, NC, NR, xc,     &
       dr_xave,                &
       rho, tem                )
    use mod_atmos_cnst, only: &
         I_QC
    implicit none
    !
    integer, intent(in)  :: ijdim
    integer, intent(in)  :: kdim
    integer, intent(in)  :: kmin
    integer, intent(in)  :: kmax 
    !
    real(8), intent(out) :: PLCaut(ijdim,kdim) ! Lc change for Auto-conversion
    real(8), intent(out) :: PNCaut(ijdim,kdim) ! Nc
    real(8), intent(out) :: PNRaut(ijdim,kdim) ! Nr 
    real(8), intent(out) :: PLCacc(ijdim,kdim) ! Lc change for Accretion
    real(8), intent(out) :: PNCacc(ijdim,kdim) ! Nc
    real(8), intent(out) :: PNRslc(ijdim,kdim) ! Nr change for Self-collection
    real(8), intent(out) :: PNRbrk(ijdim,kdim) ! Nr change for Breakup
    !
    real(8), intent(in)  :: LC(ijdim,kdim)
    real(8), intent(in)  :: LR(ijdim,kdim)
    real(8), intent(in)  :: NC(ijdim,kdim)
    real(8), intent(in)  :: NR(ijdim,kdim)
    real(8), intent(in)  :: xc(ijdim,kdim)
    real(8), intent(in)  :: dr_xave(ijdim,kdim)
    real(8), intent(in)  :: rho(ijdim,kdim)
    real(8), intent(in)  :: tem(ijdim,kdim)
    !
    ! parameter for autoconversion
    real(8), parameter :: kcc     = 4.44d9  ! collision efficiency [m3/kg2/sec]
    real(8), parameter :: tau_min = 1.d-20  ! empirical filter by T.Mitsui
    real(8), parameter :: rx_sep  = 1.d0/x_sep ! 1/x_sep, 10/08/03 [Add] T.Mitsui
    !
    ! parameter for accretion
    real(8), parameter :: kcr     = 5.8     ! collision efficiency [m3/kg2/sec]
    real(8), parameter :: thr_acc = 5.d-5   ! threshold for universal function original
    !
    ! parameter for self collection and collison break-up
    real(8), parameter :: krr     = 4.33d0  ! k_rr,      S08 (35)
    real(8), parameter :: kaprr   = 60.7d0  ! kappa_rr,  SB06(11)
    real(8), parameter :: kbr     = 1000.D0 ! k_br,      SB06(14)
    real(8), parameter :: kapbr   = 2.3d3   ! kappa_br,  SB06(15)
    real(8), parameter :: dr_min  = 0.35d-3 ! minimum diameter, SB06(13)-(15)
    !
    ! work variables
    real(8) :: coef_nuc0 ! coefficient of number for Auto-conversion
    real(8) :: coef_nuc1 !                mass   
    real(8) :: coef_aut0 !                number
    real(8) :: coef_aut1 !                mass
    real(8) :: lwc       ! lc+lr
    real(8) :: tau       ! conversion ratio: qr/(qc+qr) ranges [0:1]
    real(8) :: rho_fac   ! factor of air density 
    real(8) :: psi_aut   ! Universal function of Auto-conversion
    real(8) :: psi_acc   ! Universal function of Accretion
    real(8) :: psi_brk   ! Universal function of Breakup
    real(8) :: ddr       ! diameter difference from equilibrium 
    !
    integer :: ij, k
    !
    PLCaut(:,1:kmin)=0.D0
    PNCaut(:,1:kmin)=0.D0
    PNRaut(:,1:kmin)=0.D0
    PLCacc(:,1:kmin)=0.D0
    PNCacc(:,1:kmin)=0.D0
    PNRslc(:,1:kmin)=0.D0
    PNRbrk(:,1:kmin)=0.D0
    !
    PLCaut(:,kmax:kdim)=0.D0
    PNCaut(:,kmax:kdim)=0.D0
    PNRaut(:,kmax:kdim)=0.D0
    PLCacc(:,kmax:kdim)=0.D0
    PNCacc(:,kmax:kdim)=0.D0
    PNRslc(:,kmax:kdim)=0.D0
    PNRbrk(:,kmax:kdim)=0.D0
    !
    coef_nuc0 = (nu(I_QC)+2.d0)/(nu(I_QC)+1.d0)
    coef_nuc1 = (nu(I_QC)+2.d0)*(nu(I_QC)+4.d0)/(nu(I_QC)+1.d0)/(nu(I_QC)+1.d0)
    coef_aut0 =  -kcc*coef_nuc0
    coef_aut1 =  -kcc/x_sep/20.D0*coef_nuc1
    !
    do k=kmin, kmax
       do ij=1, ijdim
          lwc = LR(ij,k)+LC(ij,k) 
          if( lwc > xc_min )then
             tau  = max(tau_min, LR(ij,k)/lwc)
          else
             tau  = tau_min
          end if
          rho_fac = sqrt(rho_0/max(rho(ij,k),rho_min)) 
          !
          ! Auto-conversion ( cloud-cloud => rain )
          psi_aut       = 400.D0*(tau**0.7d0)*(1.d0 - (tau**0.7d0))**3   ! (6) SB06
          PNCaut(ij,k)  = coef_aut0*LC(ij,k)*LC(ij,k)*rho_fac*rho_fac    ! (9) SB06 sc+aut
          ! lc = lwc*(1-tau), lr = lwc*tau
          PLCaut(ij,k)  = coef_aut1*lwc*lwc*xc(ij,k)*xc(ij,k) &          ! (4) SB06
               *((1.d0-tau)*(1.d0-tau) + psi_aut)*rho_fac*rho_fac        !
          PNRaut(ij,k)  = -rx_sep*PLCaut(ij,k)                           ! (A7) SB01
          !
          ! Accretion ( cloud-rain => rain )
          psi_acc       =(tau/(tau+thr_acc))**4                          ! (8) SB06
          PLCacc(ij,k)  = -kcr*LC(ij,k)*LR(ij,k)*rho_fac*psi_acc         ! (7) SB06
          PNCacc(ij,k)  = -kcr*NC(ij,k)*LR(ij,k)*rho_fac*psi_acc         ! (A6) SB01
          !
          ! Self-collection ( rain-rain => rain )
          PNRslc(ij,k)  = -krr*NR(ij,k)*LR(ij,k)*rho_fac                 ! (A.8) SB01 
          !
          ! Collisional breakup of rain
          ddr           = min(1.d-3, dr_xave(ij,k) - dr_eq )
          if      (dr_xave(ij,k) < dr_min )then                          ! negligible
             psi_brk      = -1.d0
             PNRbrk(ij,k) = 0.D0
          else if (dr_xave(ij,k) <= dr_eq  )then
             psi_brk      = kbr*ddr + 1.d0                               ! (14) SB06 (+1 is necessary)
             PNRbrk(ij,k) = - (psi_brk + 1.d0)*PNRslc(ij,k)              ! (13) SB06
          else
             psi_brk      = 2.d0*exp(kapbr*ddr) - 1.d0                   ! (15) SB06
             PNRbrk(ij,k) = - (psi_brk + 1.d0)*PNRslc(ij,k)              ! (13) SB06
          end if
          !
       end do
    end do
    !
    return
  end subroutine aut_acc_slc_brk
  ! Vapor Deposition, Ice Melting
  subroutine dep_vapor_melt_ice( &
       ijdim, kdim,          & ! in
       kmin, kmax,           & ! in
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
    !
    use mod_atmos_cnst, only : &
         CNST_RVAP,      &
         CNST_RAIR,      &
         CNST_TEM00,     &
         CNST_LH0,       &
         CNST_LHS0,      &
         CNST_LHF00,     &
         CNST_PI,        &
         CNST_CP,        &
         CNST_PRES0,     &
         CNST_EMELT,     & 
         CNST_PSAT0,     &
         I_QC, I_QR, I_QI, I_QS, I_QG, &
         HYDRO_MAX, &
         LHV, LHS              
    use mod_satadjust, only : &
         moist_psat_water,    &
         moist_psat_ice
    !
    implicit none
    !
    integer, intent(in)  :: ijdim
    integer, intent(in)  :: kdim
    integer, intent(in)  :: kmin
    integer, intent(in)  :: kmax
    ! Diffusion growth or Evaporation, Sublimation
    real(8), intent(out) :: PLCdep(ijdim,kdim)  ! mass change   for cloud, [Add]  09/08/18 T.Mitsui
    real(8), intent(out) :: PLRdep(ijdim,kdim)  ! mass change   for rain deposion
    real(8), intent(out) :: PNRdep(ijdim,kdim)  ! number change 
    real(8), intent(out) :: PLIdep(ijdim,kdim)  ! mass          for cloud ice
    real(8), intent(out) :: PNIdep(ijdim,kdim)  ! number        
    real(8), intent(out) :: PLSdep(ijdim,kdim)  ! mass          for snow
    real(8), intent(out) :: PNSdep(ijdim,kdim)  ! number        
    real(8), intent(out) :: PLGdep(ijdim,kdim)  ! mass          for graupel
    real(8), intent(out) :: PNGdep(ijdim,kdim)  ! number
    ! Melting under condition(T > 273.15K and Gm > 0.0 )
    real(8), intent(out) :: PLImlt(ijdim,kdim)  ! mass          for cloud ice melting
    real(8), intent(out) :: PNImlt(ijdim,kdim)  ! number        
    real(8), intent(out) :: PLSmlt(ijdim,kdim)  ! mass          for snow
    real(8), intent(out) :: PNSmlt(ijdim,kdim)  ! number        
    real(8), intent(out) :: PLGmlt(ijdim,kdim)  ! mass          for graupel
    real(8), intent(out) :: PNGmlt(ijdim,kdim)  ! number
    !
    real(8), intent(in)  :: rho(ijdim,kdim)     ! air density
    real(8), intent(in)  :: tem(ijdim,kdim)     ! air temperature
    real(8), intent(in)  :: pre(ijdim,kdim)     ! air pressure
    real(8), intent(in)  :: qd(ijdim,kdim)      ! mixing ratio of dry air
    real(8), intent(in)  :: esw(ijdim,kdim)     ! saturation vapor pressure(liquid water)
    real(8), intent(in)  :: esi(ijdim,kdim)     ! saturation vapor pressure(solid water)
    real(8), intent(in)  :: LV(ijdim,kdim)      ! mass   of vapor
    real(8), intent(in)  :: NC(ijdim,kdim)      ! number of cloud  09/08/18 [Add] T.Mitsui
    real(8), intent(in)  :: NR(ijdim,kdim)      ! number of rain
    real(8), intent(in)  :: NI(ijdim,kdim)      !        of cloud ice
    real(8), intent(in)  :: NS(ijdim,kdim)      !        of snow
    real(8), intent(in)  :: NG(ijdim,kdim)      !        of graupel
    real(8), intent(in)  :: LI(ijdim,kdim)      ! mass   of cloud ice
    real(8), intent(in)  :: LS(ijdim,kdim)      ! mass   of snow
    real(8), intent(in)  :: LG(ijdim,kdim)      ! mass   of graupel
    real(8), intent(in)  :: xc(ijdim,kdim)      ! mean mass of cloud(filtered) [Add] 09/08/18 T.Mitsui
    real(8), intent(in)  :: xr(ijdim,kdim)      ! mean mass of rain(filtered)
    real(8), intent(in)  :: xi(ijdim,kdim)      !           of cloud ice(filtered)
    real(8), intent(in)  :: xs(ijdim,kdim)      !           of snow(filtered)
    real(8), intent(in)  :: xg(ijdim,kdim)      !           of graupel(filtered)
    ! Notice following values differ from mean terminal velocity or diameter.
    ! mean(vt(x)) /= vt(mean(x)) and mean(D(x)) /= D(mean(x)) 
    ! Following ones are vt(mean(x)) and D(mean(x)).
    real(8), intent(in)  :: vt_xave(ijdim,kdim,HYDRO_MAX,2) ! terminal velocity of mean cloud 09/08/18 [Add], T.Mitsui
    !
    real(8), intent(in)  :: dc_xave(ijdim,kdim) ! diameter of mean cloud 09/08/18 [Add] T.Mitsui
    real(8), intent(in)  :: dr_xave(ijdim,kdim) ! diameter of mean rain
    real(8), intent(in)  :: di_xave(ijdim,kdim) !                  ice
    real(8), intent(in)  :: ds_xave(ijdim,kdim) !                  snow
    real(8), intent(in)  :: dg_xave(ijdim,kdim) !                  graupel
    !
    real(8) :: rho_lim            ! limited density              09/08/18 T.Mitsui
    real(8) :: temc_lim           ! limited temperature[celsius] 09/08/18 T.Mitsui
    real(8) :: pre_lim            ! limited density              09/08/18 T.Mitsui
    real(8) :: temc               ! temperature[celsius]
    real(8) :: pv                 ! vapor pressure
    real(8) :: qv                 ! mixing ratio of water vapor [Add] 09/08/18
    real(8) :: ssw                ! super saturation ratio(liquid water)
    real(8) :: ssi                ! super saturation ratio(ice water)
    real(8) :: nua, r_nua         ! kinematic viscosity of air
    real(8) :: mua                ! viscosity of air
    real(8) :: Ka                 ! thermal conductance
    real(8) :: Dw                 ! diffusivity of water vapor
    real(8) :: Dt                 ! diffusivity of heat
    real(8) :: Gw, Gi             ! diffusion factor by balance between heat and vapor
    real(8) :: Gwr, Gii, Gis, Gig ! for rain, ice, snow and graupel.
    real(8) :: Gm                 ! melting factor by balance between heat and vapor
    real(8) :: Nsc_r3             !
    ! [Mod] 11/08/30 T.Mitsui, considering large and small branches
    real(8) :: Nrecs_r2            ! 09/08/18 [Add] T.Mitsui
    real(8) :: Nrers_r2, Nreis_r2  !
    real(8) :: Nress_r2, Nregs_r2  !
    real(8) :: Nrecl_r2            ! 09/08/18 [Add] T.Mitsui
    real(8) :: Nrerl_r2, Nreil_r2  !
    real(8) :: Nresl_r2, Nregl_r2  !
    real(8) :: NscNrer_s, NscNrer_l
    real(8) :: NscNrei_s, NscNrei_l
    real(8) :: NscNres_s, NscNres_l
    real(8) :: NscNreg_s, NscNreg_l
    real(8) :: ventLR_s, ventLR_l
    real(8) :: ventNI_s, ventNI_l, ventLI_s, ventLI_l
    real(8) :: ventNS_s, ventNS_l, ventLS_s, ventLS_l
    real(8) :: ventNG_s, ventNG_l, ventLG_s, ventLG_l
    !
    real(8) :: wtr, wti, wts, wtg
    real(8), parameter :: r_14=1.d0/1.4d0
    real(8), parameter :: r_15=1.d0/1.5d0
    !
    real(8) :: ventNR, ventLR(ijdim,kdim)     !
    real(8) :: ventNI(ijdim,kdim), ventLI(ijdim,kdim)     !
    real(8) :: ventNS(ijdim,kdim), ventLS(ijdim,kdim)     !
    real(8) :: ventNG(ijdim,kdim), ventLG(ijdim,kdim)     !
    real(8) :: ah_vent1_rs(ijdim,kdim)
    real(8) :: ah_vent1_rl(ijdim,kdim)
    real(8) :: bh_vent1_rs(ijdim,kdim)
    real(8) :: bh_vent1_rl(ijdim,kdim)
    !
    real(8), parameter :: Re_max=1.d3
    real(8), parameter :: Re_min=1.d-4
    !
    integer :: ij, k
    !
    PLCdep(:,1:kmin)=0.D0 
    PLRdep(:,1:kmin)=0.D0
    PNRdep(:,1:kmin)=0.D0
    PLIdep(:,1:kmin)=0.D0
    PNIdep(:,1:kmin)=0.D0
    PLSdep(:,1:kmin)=0.D0
    PNSdep(:,1:kmin)=0.D0
    PLGdep(:,1:kmin)=0.D0
    PNGdep(:,1:kmin)=0.D0
    PLImlt(:,1:kmin)=0.D0
    PNImlt(:,1:kmin)=0.D0
    PLSmlt(:,1:kmin)=0.D0
    PNsmlt(:,1:kmin)=0.D0
    PLGmlt(:,1:kmin)=0.D0
    PNGmlt(:,1:kmin)=0.D0
    !
    PLCdep(:,kmax:kdim)=0.D0
    PLRdep(:,kmax:kdim)=0.D0
    PNRdep(:,kmax:kdim)=0.D0
    PLIdep(:,kmax:kdim)=0.D0
    PNIdep(:,kmax:kdim)=0.D0
    PLSdep(:,kmax:kdim)=0.D0
    PNSdep(:,kmax:kdim)=0.D0
    PLGdep(:,kmax:kdim)=0.D0
    PNGdep(:,kmax:kdim)=0.D0
    PLImlt(:,kmax:kdim)=0.D0
    PNImlt(:,kmax:kdim)=0.D0
    PLSmlt(:,kmax:kdim)=0.D0
    PNsmlt(:,kmax:kdim)=0.D0
    PLGmlt(:,kmax:kdim)=0.D0
    PNGmlt(:,kmax:kdim)=0.D0
    !
    ah_vent1_rs(:,:)  = ah_vent1(I_QR,1)
    ah_vent1_rl(:,:)  = ah_vent1(I_QR,2)
    bh_vent1_rs(:,:)  = bh_vent1(I_QR,1)
    bh_vent1_rl(:,:)  = bh_vent1(I_QR,2)
    !
    ventNR=0.D0
    ventNI=0.D0
    ventNS=0.D0
    ventNG=0.D0
    ventLR=0.D0
    ventLI=0.D0
    ventLS=0.D0
    ventLG=0.D0
    !
    ! Notice,T.Mitsui
    ! Vapor deposition and melting would not be solved iteratively to reach equilibrium.
    ! Because following phenomena are not adjustment but transition.
    ! Just time-scales differ among them.
    ! If we would treat more appropreately, there would be time-splitting method to solve each ones.
    do k=kmin, kmax
       do ij=1, ijdim
          temc    = tem(ij,k) - CNST_TEM00   ! degC
          temc_lim= max(temc, -40.D0 )       ! [Add] 09/08/18 T.Mitsui, Pruppacher and Klett(1997),(13-3)
          rho_lim = max(rho(ij,k),rho_min)   ! [Add] 09/08/18 T.Mitsui
          qv      = LV(ij,k)/rho_lim
          pre_lim = rho_lim*(qd(ij,k)*CNST_RAIR + qv*CNST_RVAP)*(temc_lim+CNST_TEM00) ![Add] 09/08/18 T.Mitsui
          !--------------------------------------------------------------------
          ! Diffusion growth part is described in detail
          ! by Pruppacher and Klett (1997) Sec. 13.2(liquid) and 13.3(solid)
          !
          ! G:factor of thermal diffusion(1st.term) and vapor diffusion(2nd. term)
          ! SB06(23),(38), Lin et al(31),(52) or others
          ! Dw is introduced by Pruppacher and Klett(1997),(13-3)
          Dw      = 0.211d-4* (((temc_lim+CNST_TEM00)/CNST_TEM00)**1.94) *(CNST_PRES0/pre_lim)
          Ka      = Ka0  + temc_lim*dKa_dT  
          mua     = mua0 + temc_lim*dmua_dT 
          nua     = mua/rho_lim
          r_nua   = 1.d0/nua
          Gw      = (CNST_LH0 /Ka/tem(ij,k))*(CNST_LH0 /CNST_RVAP/tem(ij,k)-1.0D0)+(CNST_RVAP*tem(ij,k)/Dw/esw(ij,k))
          Gi      = (CNST_LHS0/Ka/tem(ij,k))*(CNST_LHS0/CNST_RVAP/tem(ij,k)-1.0D0)+(CNST_RVAP*tem(ij,k)/Dw/esi(ij,k))
          ! capacities account for their surface geometries
          Gwr     = 4.d0*CNST_PI/cap(I_QR)/Gw
          Gii     = 4.d0*CNST_PI/cap(I_QI)/Gi
          Gis     = 4.d0*CNST_PI/cap(I_QS)/Gi
          Gig     = 4.d0*CNST_PI/cap(I_QG)/Gi
          ! vent: ventilation effect( asymmetry vapor field around particles due to aerodynamic )
          ! SB06 (30),(31) and each coefficient is by (88),(89) 
          Nsc_r3  = (nua/Dw)**(0.33333333d0)                    ! (Schmidt number )^(1/3)
          !
          Nrecs_r2 = sqrt(max(Re_min,min(Re_max,vt_xave(ij,k,I_QC,1)*dc_xave(ij,k)*r_nua))) ! (Reynolds number)^(1/2) cloud
          Nrers_r2 = sqrt(max(Re_min,min(Re_max,vt_xave(ij,k,I_QR,1)*dr_xave(ij,k)*r_nua))) ! (Reynolds number)^(1/2) rain
          Nreis_r2 = sqrt(max(Re_min,min(Re_max,vt_xave(ij,k,I_QI,1)*di_xave(ij,k)*r_nua))) ! (Reynolds number)^(1/2) cloud ice
          Nress_r2 = sqrt(max(Re_min,min(Re_max,vt_xave(ij,k,I_QS,1)*ds_xave(ij,k)*r_nua))) ! (Reynolds number)^(1/2) snow
          Nregs_r2 = sqrt(max(Re_min,min(Re_max,vt_xave(ij,k,I_QG,1)*dg_xave(ij,k)*r_nua))) ! (Reynolds number)^(1/2) graupel    
          !
          Nrecl_r2 = sqrt(max(Re_min,min(Re_max,vt_xave(ij,k,I_QC,2)*dc_xave(ij,k)*r_nua))) ! (Reynolds number)^(1/2) cloud
          Nrerl_r2 = sqrt(max(Re_min,min(Re_max,vt_xave(ij,k,I_QR,2)*dr_xave(ij,k)*r_nua))) ! (Reynolds number)^(1/2) rain
          Nreil_r2 = sqrt(max(Re_min,min(Re_max,vt_xave(ij,k,I_QI,2)*di_xave(ij,k)*r_nua))) ! (Reynolds number)^(1/2) cloud ice
          Nresl_r2 = sqrt(max(Re_min,min(Re_max,vt_xave(ij,k,I_QS,2)*ds_xave(ij,k)*r_nua))) ! (Reynolds number)^(1/2) snow
          Nregl_r2 = sqrt(max(Re_min,min(Re_max,vt_xave(ij,k,I_QG,2)*dg_xave(ij,k)*r_nua))) ! (Reynolds number)^(1/2) graupel    
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
          wtr     = ( min(max( NscNrer_s*r_14, 0.5d0), 2.0d0) -0.5d0 )*r_15 ! weighting between 1.4*0.5 and 1.4*2
          wti     = ( min(max( NscNrei_s     , 0.5d0), 2.0d0) -0.5d0 )*r_15 ! weighting between 1.0*0.5 and 1.0*2
          wts     = ( min(max( NscNres_s*r_14, 0.5d0), 2.0d0) -0.5d0 )*r_15 ! weighting between 1.4*0.5 and 1.4*2
          wtg     = ( min(max( NscNreg_s*r_14, 0.5d0), 2.0d0) -0.5d0 )*r_15 ! weighting between 1.4*0.5 and 1.4*2
          ! interpolation between two branches
          ventNI(ij,k)  = (1.d0-wti)*ventNI_s + wti*ventNI_l
          ventNS(ij,k)  = (1.d0-wts)*ventNS_s + wts*ventNS_l
          ventNG(ij,k)  = (1.d0-wtg)*ventNG_s + wtg*ventNG_l
          !
          ventLR(ij,k)  = (1.d0-wtr)*ventLR_s + wtr*ventLR_l
          ventLI(ij,k)  = (1.d0-wti)*ventLI_s + wti*ventLI_l
          ventLS(ij,k)  = (1.d0-wts)*ventLS_s + wts*ventLS_l
          ventLG(ij,k)  = (1.d0-wtg)*ventLG_s + wtg*ventLG_l
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
!!$       PLCdep(ij,k) = Gwr*ssw*NC(ij,k)*dc_xave(ij,k)*coef_deplc 
!!$       PLRdep(ij,k) = Gwr*ssw*NR(ij,k)*dr_xave(ij,k)*ventLR     
!!$       PLIdep(ij,k) = Gii*ssi*NI(ij,k)*di_xave(ij,k)*ventLI
!!$       PLSdep(ij,k) = Gis*ssi*NS(ij,k)*ds_xave(ij,k)*ventLS
!!$       PLGdep(ij,k) = Gig*ssi*NG(ij,k)*dg_xave(ij,k)*ventLG
          PLCdep(ij,k) = Gwr*NC(ij,k)*dc_xave(ij,k)*coef_deplc 
          PLRdep(ij,k) = Gwr*NR(ij,k)*dr_xave(ij,k)*ventLR(ij,k)     
          PLIdep(ij,k) = Gii*NI(ij,k)*di_xave(ij,k)*ventLI(ij,k)
          PLSdep(ij,k) = Gis*NS(ij,k)*ds_xave(ij,k)*ventLS(ij,k)
          PLGdep(ij,k) = Gig*NG(ij,k)*dg_xave(ij,k)*ventLG(ij,k)
          PNRdep(ij,k) = PLRdep(ij,k)/xr(ij,k) 
          PNIdep(ij,k) = 0.D0
          PNSdep(ij,k) = PLSdep(ij,k)/xs(ij,k) 
          PNGdep(ij,k) = PLGdep(ij,k)/xg(ij,k) 
          !
          !------------------------------------------------------------------------
          ! Melting part is described by Pruppacher and Klett (1997) Sec.16.3.1
          ! Here we omit "Shedding" of snow-flakes and ice-particles.
          ! "Shedding" may be applicative if you refer 
          ! eq.(38) in Cotton etal.(1986) Jour. Clim. Appl. Meteor. p.1658-1680.
          ! SB06(73)
          Dt      = Ka/(CNST_CP*rho_0)
          ! Gm: factor caused by balance between 
          !     "water evaporation cooling(1st.)" and "fusion heating(2nd.)" 
          ! SB06(76)
          ! [fix] 08/05/08 T.Mitsui  LHF00 => EMELT  and  esw => PSAT0
          ! CNST_LHS0 is more suitable than LHS because melting occurs around 273.15 K.
          Gm      = 2.d0*CNST_PI/CNST_EMELT&
               * ( (Ka*Dt/Dw)*(temc) + (Dw*CNST_LHS0/CNST_RVAP)*(esi(ij,k)/tem(ij,k)-CNST_PSAT0/CNST_TEM00) )
          ! SB06(76)
          ! Notice! melting only occurs where T > 273.15 K else doesn't.
          ! [fix] 08/05/08 T.Mitsui, Gm could be both positive and negative value.
          !       See Pruppacher and Klett(1997) eq.(16-79) or Rasmussen and Pruppacher(1982)
          if( (temc>=0.D0) .and. (Gm>0.D0) )then !  if Gm==0 then rh and tem is critical value for melting process.
             ! 08/05/16 [Mod] T.Mitsui, change term of PLimlt. N_i => L_i/ (limited x_i)
             ! because melting never occur when N_i=0.
             PLImlt(ij,k) = - Gm * LI(ij,k)*di_xave(ij,k)*ventLI(ij,k)/xi(ij,k)
             ! [Mod] 08/08/23 T.Mitsui for Seifert(2008)
             PNImlt(ij,k) = - Gm * NI(ij,k)*di_xave(ij,k)*ventNI(ij,k)/xi(ij,k) ! 09/08/18 [Mod] recover, T.Mitsui
             PLSmlt(ij,k) = - Gm * LS(ij,k)*ds_xave(ij,k)*ventLS(ij,k)/xs(ij,k)
             ! [Mod] 08/08/23 T.Mitsui for Seifert(2008)
             PNSmlt(ij,k) = - Gm * NS(ij,k)*ds_xave(ij,k)*ventNS(ij,k)/xs(ij,k) ! 09/08/18 [Mod] recover, T.Mitsui
             PLGmlt(ij,k) = - Gm * LG(ij,k)*dg_xave(ij,k)*ventLG(ij,k)/xg(ij,k)
             ! [Mod] 08/08/23 T.Mitsui for Seifert(2008)
             PNGmlt(ij,k) = - Gm * NG(ij,k)*dg_xave(ij,k)*ventNG(ij,k)/xg(ij,k) ! 09/08/18 [Mod] recover, T.Mitsui
          else
             PLImlt(ij,k) = 0.D0
             PNImlt(ij,k) = 0.D0
             PLSmlt(ij,k) = 0.D0
             PNSmlt(ij,k) = 0.D0
             PLGmlt(ij,k) = 0.D0
             PNGmlt(ij,k) = 0.D0
          end if
          !
       end do
    end do
    !
    return
  end subroutine dep_vapor_melt_ice
  !
  subroutine freezing_water( &
       ijdim, kdim,          &
       kmin, kmax,           &
       dt,                   & 
       PLChom, PNChom,       &
       PLChet, PNChet,       &
       PLRhet, PNRhet,       &
       LC, LR, NC, NR, xc, xr, tem   )
    !
    ! In this subroutine, 
    ! We assumed surface temperature of droplets are same as environment.
    use mod_atmos_cnst, only : &
         CNST_TEM00,     &
         I_QC, I_QR
    implicit none
    !
    integer, intent(in) :: ijdim
    integer, intent(in) :: kdim
    integer, intent(in) :: kmin
    integer, intent(in) :: kmax
    !
    real(8), intent(in) :: dt  
    ! freezing in nature (homogenous)
    real(8), intent(out):: PLChom(ijdim,kdim) ! cloud water => cloud ice
    real(8), intent(out):: PNChom(ijdim,kdim) ! cloud water => cloud ice
    ! freezing via aerosols (heterogenous)
    real(8), intent(out):: PLChet(ijdim,kdim) ! cloud water => cloud ice
    real(8), intent(out):: PLRhet(ijdim,kdim) ! rain        => graupel
    real(8), intent(out):: PNChet(ijdim,kdim) ! cloud water => cloud ice
    real(8), intent(out):: PNRhet(ijdim,kdim) ! rain        => graupel
    !
    real(8), intent(in) :: tem(ijdim,kdim)
    !
    real(8), intent(in) :: LC(ijdim,kdim)
    real(8), intent(in) :: LR(ijdim,kdim)
    real(8), intent(in) :: NC(ijdim,kdim) 
    real(8), intent(in) :: NR(ijdim,kdim) 
    real(8), intent(in) :: xc(ijdim,kdim)
    real(8), intent(in) :: xr(ijdim,kdim)
    !
    real(8), parameter :: temc_min = -65.d0
    real(8), parameter :: a_het = 0.2d0  ! SB06 (44)
    real(8), parameter :: b_het = 0.65d0 ! SB06 (44)
    !
    real(8) :: coef_m2_c
    real(8) :: coef_m2_r
    ! temperature [celsius]
    real(8) :: temc, temc2, temc3, temc4
    ! temperature function of homegenous/heterogenous freezing
    real(8) :: Jhom, Jhet
    real(8) :: rdt
    !
    integer :: ij,k
    !
    PLChom(:,1:kmin)=0.D0
    PNChom(:,1:kmin)=0.D0
    PLChet(:,1:kmin)=0.D0
    PNChet(:,1:kmin)=0.D0
    PLRhet(:,1:kmin)=0.D0
    PNRhet(:,1:kmin)=0.D0
    !
    PLChom(:,kmax:kdim)=0.D0
    PNChom(:,kmax:kdim)=0.D0
    PLChet(:,kmax:kdim)=0.D0
    PNChet(:,kmax:kdim)=0.D0
    PLRhet(:,kmax:kdim)=0.D0
    PNRhet(:,kmax:kdim)=0.D0
    !
    rdt = 1.d0/dt
    !
    coef_m2_c =   coef_m2(I_QC)
    coef_m2_r =   coef_m2(I_QR)
    !
    do k=kmin, kmax
       do ij=1, ijdim
          temc = max( tem(ij,k) - CNST_TEM00, temc_min )
          ! These cause from aerosol-droplet interaction.
          ! Bigg(1953) formula, Khain etal.(2000) eq.(4.5), Pruppacher and Klett(1997) eq.(9-48)
          Jhet =  a_het*exp( -b_het*temc - 1.d0 )
          ! These cause in nature.
          ! Cotton and Field 2002, QJRMS. (12)
          if( temc < -65.d0 )then
             jhom = 10.D0**(24.37236d0)*1.d3
             Jhet =  a_het*exp( 65.d0*b_het - 1.d0 ) ! 09/04/14 [Add], fixer T.Mitsui
          else if( temc < -30.D0 ) then
             temc2 = temc*temc
             temc3 = temc*temc2
             temc4 = temc2*temc2
             jhom = 10.D0**(&
                  - 243.4d0 - 14.75d0*temc - 0.307d0*temc2 &
                  - 0.00287d0*temc3 - 0.0000102*temc4 ) *1.d3
          else if( temc < 0.D0) then
             jhom = 10.D0**(-7.63d0-2.996d0*(temc+30.D0))*1.d3
          else
             Jhom = 0.D0
             Jhet = 0.D0
          end if
          ! Note, xc should be limited in range[xc_min:xc_max].
          ! and PNChom need to be calculated by NC
          ! because reduction rate of Nc need to be bound by NC.
          ! For the same reason PLChom also bound by LC and xc.
          ! Basically L and N should be independent 
          ! but average particle mass x should be in suitable range.
          ! Homogenous Freezing
          PLChom(ij,k) = 0.D0
          PNChom(ij,k) = 0.D0
          ! Heterogenous Freezing
          PLChet(ij,k) = -rdt*LC(ij,k)*( 1.d0 - exp( -coef_m2_c*xc(ij,k)*(Jhet+Jhom)*dt ) )
          PNChet(ij,k) = -rdt*NC(ij,k)*( 1.d0 - exp( -          xc(ij,k)*(Jhet+Jhom)*dt ) )    
          PLRhet(ij,k) = -rdt*LR(ij,k)*( 1.d0 - exp( -coef_m2_r*xr(ij,k)*(Jhet+Jhom)*dt ) )
          PNRhet(ij,k) = -rdt*NR(ij,k)*( 1.d0 - exp( -          xr(ij,k)*(Jhet+Jhom)*dt ) )    
       end do
    end do
    !
    return
  end subroutine freezing_water
  !
  function gammafunc( xx ) result(f)
    implicit none
    real(8), intent(in) :: xx
    real(8) :: f
    real(8) :: coef(6)=(/&
         +76.18009172947146D0,&
         -86.50532032941677D0,&
         +24.01409824083091D0,&
         -1.231739572450155D0,&
         +0.1208650973866179D-2,&
         -0.5395239384953D-5&
         /)
    integer :: j
    real(8) :: x,y,tmp,ser
    
    x=xx
    y=x
    tmp=x+5.5D0
    tmp = tmp - (x+0.5)*log(tmp)
    ser=1.000000000190015D0
    do j=1,6
       y=y+1
       ser = ser+coef(j)/y
    end do
    f = exp(-tmp+log(2.5066282746310005D0*ser/x))
  end function gammafunc
  !
  function gammafunc_3d( ijdim, kdim, x ) 
    implicit none
    integer, intent(in) :: ijdim
    integer, intent(in) :: kdim
    real(8), intent(in) :: x(ijdim,kdim)
    real(8)             :: gammafunc_3d(ijdim,kdim)
    real(8), parameter  :: coef(6)=(/&
         +76.18009172947146D0,&
         -86.50532032941677D0,&
         +24.01409824083091D0,&
         -1.231739572450155D0,&
         +0.1208650973866179D-2,&
         -0.5395239384953D-5&
         /)
    real(8), parameter :: ser0=1.000000000190015D0
    real(8) :: tmp(ijdim,kdim)
    real(8) :: ser(ijdim,kdim)
    integer :: ij,k,iter
    !
    ser(:,:) = ser0 &
         + coef(1)/(x(:,:)+1.d0) + coef(2)/(x(:,:)+2.d0) + coef(3)/(x(:,:)+3.d0) &
         + coef(4)/(x(:,:)+4.d0) + coef(5)/(x(:,:)+5.d0) + coef(6)/(x(:,:)+6.d0)
    tmp(:,:) = x(:,:)+5.5d0 - (x(:,:)+0.5d0)*log(x(:,:)+5.5d0)
    gammafunc_3d(:,:)   = exp(-tmp(:,:)+log(2.5066282746310005D0*ser(:,:)/x(:,:)))
    return
  end function gammafunc_3d
  !-----------------------------------------------------------------------------
  function igammafunc_3d( ijdim, kdim, x, alpha, gm ) result(igm)
    ! Section 6.2 in "Numerical Recipes in C"
    ! function is represented by(6.2.1)
    ! g(x,alpha)=1/gamma(alpha)*integral_0^x { exp(-t) * t^(alpha-1) }dt 
    ! g(0)=0 and g(+infinity)=1
    implicit none
    !
    integer, intent(in) :: ijdim
    integer, intent(in) :: kdim              
    real(8), intent(in) :: x(ijdim,kdim)     !
    real(8), intent(in) :: alpha(ijdim,kdim) ! 
    real(8), intent(in) :: gm(ijdim,kdim)    ! gamma function
    ! incomplete gamma function
    real(8) :: igm(ijdim,kdim)
    ! 
    real(8) :: lx(ijdim,kdim)
    real(8) :: lgm(ijdim,kdim)
    ! work
    ! coefficient of expansion using in calculation of igm
    real(8) :: a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15
    real(8) :: an1,an2,an3,an4,an5
    real(8) :: b0,b1,b2,b3,b4,b5
    real(8) :: c0,c1,c2,c3,c4,c5
    real(8) :: d0,d1,d2,d3,d4,d5
    real(8) :: e0,e1,e2,e3,e4,e5
    real(8) :: h0,h1,h2,h3,h4,h5
    real(8), parameter :: eps=1.d-30
    !
    integer :: ij,k
    !
    lgm(:,:)=log(gm(:,:))
    lx(:,:) =log( x(:,:))
    !
    ! Incomplete Gamma Function
    !
    do k=1,kdim
       do ij=1,ijdim
          if     ( x(ij,k) < 1.d-2*alpha(ij,k) )then ! negligible 
             igm(ij,k)=0.D0
          else if( x(ij,k) < alpha(ij,k)+1.d0 )then ! series expansion (6.2.5) 
             !
             ! 10th-truncation is enough for cloud droplet.
             a0   = 1.d0/alpha(ij,k)         ! n=0
             a1   = a0*x(ij,k)/(alpha(ij,k)+1.d0)  ! n=1
             a2   = a1*x(ij,k)/(alpha(ij,k)+2.d0)  ! n=2
             a3   = a2*x(ij,k)/(alpha(ij,k)+3.d0)  ! n=3
             a4   = a3*x(ij,k)/(alpha(ij,k)+4.d0)  ! n=4
             a5   = a4*x(ij,k)/(alpha(ij,k)+5.d0)  ! n=5
             !
             a6   = a5*x(ij,k)/(alpha(ij,k)+6.d0)  ! n=6
             a7   = a6*x(ij,k)/(alpha(ij,k)+7.d0)  ! n=7
             a8   = a7*x(ij,k)/(alpha(ij,k)+8.d0)  ! n=8
             a9   = a8*x(ij,k)/(alpha(ij,k)+9.d0)  ! n=9
             a10  = a9*x(ij,k)/(alpha(ij,k)+10.D0) ! n=10
             !
             a11  = a10*x(ij,k)/(alpha(ij,k)+11.d0) ! n=11
             a12  = a11*x(ij,k)/(alpha(ij,k)+12.d0) ! n=12
             a13  = a12*x(ij,k)/(alpha(ij,k)+13.d0) ! n=13
             a14  = a13*x(ij,k)/(alpha(ij,k)+14.d0) ! n=14
             a15  = a14*x(ij,k)/(alpha(ij,k)+15.d0) ! n=15
             !
             igm(ij,k) = (a0+a1+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11+a12+a13+a14+a15)&
                  * exp( -x(ij,k) + alpha(ij,k)*lx(ij,k) - lgm(ij,k) )
          else if( x(ij,k) < alpha(ij,k)*1.d2 ) then ! continued fraction expansion (6.2.6)
             ! 
             ! 3rd-truncation is enough for rain droplet.
             ! setup
             b0   = x(ij,k)+1.d0-alpha(ij,k)
             c0   = 1.d0/eps
             d0   = 1.d0/b0
             h0   = d0
             ! n=1
             an1  = -     (1.d0-alpha(ij,k))
             b1   = b0 + 2.d0
             d1   = 1.d0/(an1*d0+b1)
             c1   = b1+an1/c0
             e1   = d1*c1
             h1   = h0*e1
             ! n=2
             an2  = -2.d0*(2.d0-alpha(ij,k))
             b2   = b1 + 2.d0
             d2   = 1.d0/(an2*d1+b2)
             c2   = b2+an2/c1
             e2   = d2*c2
             h2   = h1*e2
             ! n=3
             an3  = -3.d0*(3.d0-alpha(ij,k))
             b3   = b2 + 2.d0
             d3   = 1.d0/(an3*d2+b3)
             c3   = b3+an3/c2
             e3   = d3*c3
             h3   = h2*e3
             ! n=4
             an4  = -4.d0*(4.d0-alpha(ij,k))
             b4   = b3 + 2.d0
             d4   = 1.d0/(an4*d3+b4)
             c4   = b4+an4/c3
             e4   = d4*c4
             h4   = h3*e4
             ! n=5
             an5  = -5.d0*(5.d0-alpha(ij,k))
             b5   = b4 + 2.d0
             d5   = 1.d0/(an5*d4+b5)
             c5   = b5+an5/c4
             e5   = d5*c5
             h5   = h4*e5             
             !             
             igm(ij,k)  = 1.d0 - exp( -x(ij,k) + alpha(ij,k)*lx(ij,k) - lgm(ij,k) )*h5
          else   ! negligible
             igm(ij,k)  = 1.d0
          end if
       end do
    end do
    !
    return
  end function igammafunc_3d
  !-------------------------------------------------------------------------------
  function betafunc_3d( ijdim, kdim, x, w ) 
    implicit none
    integer, intent(in) :: ijdim
    integer, intent(in) :: kdim
    real(8), intent(in) :: x(ijdim,kdim)
    real(8), intent(in) :: w(ijdim,kdim)
    real(8)             :: betafunc_3d(ijdim,kdim)
    real(8), parameter  :: coef(6)=(/&
         +76.18009172947146D0,&
         -86.50532032941677D0,&
         +24.01409824083091D0,&
         -1.231739572450155D0,&
         +0.1208650973866179D-2,&
         -0.5395239384953D-5&
         /)
    real(8), parameter :: ser0=1.000000000190015D0
    real(8) :: y(ijdim,kdim)
    real(8) :: tmp_x(ijdim,kdim)
    real(8) :: tmp_w(ijdim,kdim)
    real(8) :: tmp_xw(ijdim,kdim)
    real(8) :: ser_x(ijdim,kdim)
    real(8) :: ser_w(ijdim,kdim)
    real(8) :: ser_xw(ijdim,kdim)
    real(8) :: lg_x(ijdim,kdim)
    real(8) :: lg_w(ijdim,kdim)
    real(8) :: lg_xw(ijdim,kdim)
    integer :: ij,k,iter
    !
    ! log(gamma(x))
    !
    ser_x(:,:) = ser0 &
         + coef(1)/(x(:,:)+1.d0) + coef(2)/(x(:,:)+2.d0) + coef(3)/(x(:,:)+3.d0) &
         + coef(4)/(x(:,:)+4.d0) + coef(5)/(x(:,:)+5.d0) + coef(6)/(x(:,:)+6.d0)
    tmp_x(:,:) = x(:,:)+5.5d0 - (x(:,:)+0.5d0)*log(x(:,:)+5.5d0)
    lg_x(:,:)  = -tmp_x(:,:)+log(2.5066282746310005D0*ser_x(:,:)/x(:,:))
    !
    ! log(gamma(w))
    !
    ser_w(:,:) = ser0 &
         + coef(1)/(w(:,:)+1.d0) + coef(2)/(w(:,:)+2.d0) + coef(3)/(w(:,:)+3.d0) &
         + coef(4)/(w(:,:)+4.d0) + coef(5)/(w(:,:)+5.d0) + coef(6)/(w(:,:)+6.d0)
    tmp_w(:,:) = w(:,:)+5.5d0 - (w(:,:)+0.5d0)*log(w(:,:)+5.5d0)
    lg_w(:,:)  = -tmp_w(:,:)+log(2.5066282746310005D0*ser_w(:,:)/w(:,:))
    !
    ! log(gamma(x+w))
    !
    y(:,:) = x(:,:) + w(:,:)
    ser_xw(:,:) = ser0 &
         + coef(1)/(y(:,:)+1.d0) + coef(2)/(y(:,:)+2.d0) + coef(3)/(y(:,:)+3.d0) &
         + coef(4)/(y(:,:)+4.d0) + coef(5)/(y(:,:)+5.d0) + coef(6)/(y(:,:)+6.d0)
    tmp_xw(:,:) = y(:,:)+5.5d0 - (y(:,:)+0.5d0)*log(y(:,:)+5.5d0)
    lg_xw(:,:)  = -tmp_xw(:,:)+log(2.5066282746310005D0*ser_xw(:,:)/y(:,:))
    !
    betafunc_3d(:,:) = exp( lg_x(:,:) + lg_w(:,:) - lg_xw(:,:))
    !
    return
  end function betafunc_3d

  !-----------------------------------------------------------------------------
  subroutine terminal_velocity( &
      velw, &
      rhoq, &
      temp, &
      pres  )
    use mod_const, only: &
       PI => CONST_PI
    use mod_grid, only: &
       KA   => GRID_KA,   &
       IA   => GRID_IA,   &
       JA   => GRID_JA,   &
       KS   => GRID_KS,   &
       KE   => GRID_KE,   &
       IS   => GRID_IS,   &
       IE   => GRID_IE,   &
       JS   => GRID_JS,   &
       JE   => GRID_JE
    use mod_atmos_vars, only: &
       var => atmos_var, &
       QA => A_QA, &
       I_DENS, &
       I_QC, &
       I_QR, &
       I_QI, &
       I_QS, &
       I_QG, &
       I_NC, &
       I_NR, &
       I_NI, &
       I_NS, &
       I_NG
    implicit none

    real(8), intent(out) :: velw(KA,IA,JA,QA) ! terminal velocity of cloud mass
    real(8), intent(in)  :: rhoq(KA,IA,JA,QA) ! rho * q
    real(8), intent(in)  :: temp(KA,IA,JA)    ! temperature
    real(8), intent(in)  :: pres(KA,IA,JA)    ! pressure

    real(8) :: xq      (KA,6)  ! average mass of 1 particle( mass/number )

    real(8) :: rhofac  (KA)    ! density factor for terminal velocity( air friction )
    real(8) :: rhofac_q(KA,6)

    real(8) :: rlambdar(KA)    ! work for diagnosis of Rain DSD ( Seifert, 2008 )
    real(8) :: mud_r           !
    real(8) :: dq      (KA,QA) ! work for Rogers etal. (1993)
    real(8) :: weight  (KA,QA) !
    real(8) :: velq_s  (KA,QA) ! terminal velocity for small branch of Rogers formula
    real(8) :: velq_l  (KA,QA) ! terminal velocity for large branch of Rogers formula

    integer :: k, i, j
    !---------------------------------------------------------------------------

    mud_r = 3.D0 * nu(I_QR) + 2.D0

    !OCL NORECURRENCE
    do j = JS, JE
    do i = IS, IE

       do k = KS, KE
          rhofac(k) = rho_0 / max( var(k,i,j,I_DENS), rho_min )
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
                      * ( (mud_r+3.D0) * (mud_r+2.D0) * (mud_r+1.D0) )**(-0.333333333D0)
       enddo

       ! Improved Rogers formula by T.Mitsui
       ! weigthed diameter
       do k = KS, KE
          dq(k,I_QR) = ( 4.D0 + mud_r ) * rlambdar(k) ! D^(3)+mu weighted mean diameter
          dq(k,I_NR) = ( 1.D0 + mud_r ) * rlambdar(k)
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
          velq_s(k,I_QR) = coef_vtr_ar2 * dq(k,I_QR) * ( 1.D0 - ( 1.D0 + coef_vtr_br2*rlambdar(k) )**(-5-mud_r) )
          velq_l(k,I_QR) = coef_vtr_ar1 - coef_vtr_br1 * ( 1.D0 + coef_vtr_cr1*rlambdar(k) )**(-4-mud_r)
          velq_s(k,I_NR) = coef_vtr_ar2 * dq(k,I_QR) * ( 1.D0 - ( 1.D0 + coef_vtr_br2*rlambdar(k) )**(-2-mud_r) )
          velq_l(k,I_NR) = coef_vtr_ar1 - coef_vtr_br1 * ( 1.D0 + coef_vtr_cr1*rlambdar(k) )**(-1-mud_r)
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
          weight(k,I_QR) = 0.5D0 * ( 1.D0 + tanh( PI * log( dq(k,I_QR)/d_vtr_branch ) ) ) ! Lr
          weight(k,I_QI) = 0.5D0 * ( 1.D0 + log( dq(k,I_QI)/d0_li ) )
          weight(k,I_QS) = 0.5D0 * ( 1.D0 + log( dq(k,I_QS)/d0_ls ) )
          weight(k,I_QG) = 0.5D0 * ( 1.D0 + log( dq(k,I_QG)/d0_lg ) )
          weight(k,I_NR) = 0.5D0 * ( 1.D0 + tanh( PI * log( dq(k,I_NR)/d_vtr_branch ) ) ) ! Nr
          weight(k,I_NI) = 0.5D0 * ( 1.D0 + log( dq(k,I_NI)/d0_ni ) )
          weight(k,I_NS) = 0.5D0 * ( 1.D0 + log( dq(k,I_NS)/d0_ns ) )
          weight(k,I_NG) = 0.5D0 * ( 1.D0 + log( dq(k,I_NG)/d0_ng ) )
       enddo

          ! filter is used to avoid unrealistic terminal velocity( when ni,ns,ng are too small )
       do k = KS, KE
          weight(k,I_QI) = min( max( weight(k,I_QI), 0.D0 ), 1.D0 )
          weight(k,I_QS) = min( max( weight(k,I_QS), 0.D0 ), 1.D0 )
          weight(k,I_QG) = min( max( weight(k,I_QG), 0.D0 ), 1.D0 )
          weight(k,I_NI) = min( max( weight(k,I_NI), 0.D0 ), 1.D0 )
          weight(k,I_NS) = min( max( weight(k,I_NS), 0.D0 ), 1.D0 )
          weight(k,I_NG) = min( max( weight(k,I_NG), 0.D0 ), 1.D0 )
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
       do k = KS, KE
          velw(k,i,j,I_QC) = -rhofac_q(k,I_QC) * coef_vt1(I_QC,1) * xq(k,I_QC)**beta_v(I_QC,1)
          velw(k,i,j,I_NC) = -rhofac_q(k,I_QC) * coef_vt0(I_QC,1) * xq(k,I_QC)**beta_v(I_QC,1)
       enddo
       do k = KS, KE
          velw(k,i,j,I_QR) = -rhofac_q(k,I_QR) * ( velq_l(k,I_QR) * (        weight(k,I_QR) ) &
                                                 + velq_s(k,I_QR) * ( 1.D0 - weight(k,I_QR) ) )
          velw(k,i,j,I_NR) = -rhofac_q(k,I_QR) * ( velq_l(k,I_NR) * (        weight(k,I_NR) ) &
                                                 + velq_s(k,I_NR) * ( 1.D0 - weight(k,I_NR) ) )
       enddo
       do k = KS, KE
          velw(k,i,j,I_QI) = -rhofac_q(k,I_QI) * ( velq_l(k,I_QI) * (        weight(k,I_QI) ) &
                                                 + velq_s(k,I_QI) * ( 1.D0 - weight(k,I_QI) ) )
          velw(k,i,j,I_NI) = -rhofac_q(k,I_QI) * ( velq_l(k,I_NI) * (        weight(k,I_NI) ) &
                                                 + velq_s(k,I_NI) * ( 1.D0 - weight(k,I_NI) ) )
       enddo
       do k = KS, KE
          velw(k,i,j,I_QS) = -rhofac_q(k,I_QS) * ( velq_l(k,I_QS) * (        weight(k,I_QS) ) &
                                                 + velq_s(k,I_QS) * ( 1.D0 - weight(k,I_QS) ) )
          velw(k,i,j,I_NS) = -rhofac_q(k,I_QS) * ( velq_l(k,I_NS) * (        weight(k,I_NS) ) &
                                                 + velq_s(k,I_NS) * ( 1.D0 - weight(k,I_NS) ) )
       enddo
       do k = KS, KE
          velw(k,i,j,I_QG) = -rhofac_q(k,I_QG) * ( velq_l(k,I_QG) * (        weight(k,I_QG) ) &
                                                 + velq_s(k,I_QG) * ( 1.D0 - weight(k,I_QG) ) )
          velw(k,i,j,I_NG) = -rhofac_q(k,I_QG) * ( velq_l(k,I_NG) * (        weight(k,I_NG) ) &
                                                 + velq_s(k,I_NG) * ( 1.D0 - weight(k,I_NG) ) )
       enddo

    enddo
    enddo

    return
  end subroutine terminal_velocity

  subroutine update_by_phase_change(   &
       ntdiv    , ntmax,         & ! in [Add] 10/08/03
       ijdim, kdim, kmin, kmax,  & ! in
       nqmax,                    & ! in 
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
    !
    use mod_stdio, only: &
       IO_FID_CONF,  &
       IO_FID_LOG,   &
       IO_L
    use mod_atmos_cnst, only : &
         LHV,LHF,           & 
         CVW,               &
         NNW_STR, NNW_END,  &
         NQW_STR, NQW_END,  &
         I_QV, I_QC, I_QR, I_QI, I_QS, I_QG,  &
         I_NC, I_NR, I_NI, I_NS, I_NG, &
         CNST_LH00, &
         CNST_LHS00,&
         CNST_LHF00,&
         CNST_CVV,  &
         CNST_CV,   & ! [add] 11/08/30 T.Mitsui
         CNST_CL,   &
         CNST_CI,   &
         CNST_RVAP, &
         CNST_RAIR, &
         CNST_TEM00,&
         CNST_EGRAV
    use mod_thrmdyn, only: &
         thrmdyn_qd, &
         thrmdyn_cv, &
         thrmdyn_cp
    use mod_satadjust, only : &
         moist_qsat_water, &
         moist_qsat_ice  , &
         moist_dqsw_dtem_rho, &
         moist_dqsi_dtem_rho, &
         moist_dqsw_dtem_dpre, &
         moist_dqsi_dtem_dpre
!!$    use mod_history, only: &
!!$         history_in
    !
    implicit none
    !
    integer, intent(in)    :: ntdiv                   ! [Add] 10/08/03
    integer, intent(in)    :: ntmax                   ! [Add] 10/08/03
    !
    integer, intent(in)    :: ijdim                   !
    integer, intent(in)    :: kdim, kmin, kmax        !
    integer, intent(in)    :: nqmax                   ! tracer number
    real(8), intent(in)    :: dt                      ! time step[s]
    real(8), intent(in)    :: gsgam2(ijdim,kdim)      ! metric
    real(8), intent(in)    :: z(kdim)           ! altitude [m]
    real(8), intent(in)    :: dz(kdim)          ! altitude [m]
    real(8), intent(in)    :: wh(ijdim,kdim)          ! vertical velocity @ half point[m/s]
    real(8), intent(in)    :: dTdt_rad(ijdim,kdim)    ! temperture tendency by radiation[K/s]
    real(8), intent(in)    :: rhog(ijdim,kdim)        ! density[kg/m3]
    real(8), intent(inout) :: rhoge(ijdim,kdim)       ! internal energy[J/m3]
    real(8), intent(inout) :: rhogq(ijdim,kdim,nqmax) ! tracers[kg/m3]
    real(8), intent(inout) :: q(ijdim,kdim,nqmax)     ! tracers mixing ratio[kg/kg]
    real(8), intent(inout) :: tem(ijdim,kdim)         ! temperature[K]
    real(8), intent(inout) :: pre(ijdim,kdim)         ! pressure[Pa]
    real(8), intent(in)    :: rho(ijdim,kdim)         ! air density[kg/m3]
    real(8), intent(out)   :: cva(ijdim,kdim)         ! specific heat at constant volume
    real(8), intent(in)    :: esw(ijdim,kdim)         ! saturated vapor pressure for liquid
    real(8), intent(in)    :: esi(ijdim,kdim)         !                          for ice
    real(8), intent(in)    :: lv(ijdim,kdim)          ! vapor mass [kg/m3]
    real(8), intent(in)    :: lc(ijdim,kdim), nc(ijdim,kdim) ! cloud mass [kg/m3], number[/m3]
    real(8), intent(in)    :: lr(ijdim,kdim), nr(ijdim,kdim) ! rain
    real(8), intent(in)    :: li(ijdim,kdim), ni(ijdim,kdim) ! ice
    real(8), intent(in)    :: ls(ijdim,kdim), ns(ijdim,kdim) ! snow
    real(8), intent(in)    :: lg(ijdim,kdim), ng(ijdim,kdim) ! graupel
    !+++ Freezing tendency[kg/m3/s]
    real(8), intent(inout) :: PLChom(ijdim,kdim), PNChom(ijdim,kdim) 
    real(8), intent(inout) :: PLChet(ijdim,kdim), PNChet(ijdim,kdim) 
    real(8), intent(inout) :: PLRhet(ijdim,kdim), PNRhet(ijdim,kdim) 
    !+++ Condensation/Evaporation, Deposition/Sublimaion tendency[kg/m3/s]
    real(8), intent(inout) :: PLCdep(ijdim,kdim) 
    real(8), intent(inout) :: PLRdep(ijdim,kdim), PNRdep(ijdim,kdim)
    real(8), intent(inout) :: PLIdep(ijdim,kdim), PNIdep(ijdim,kdim)
    real(8), intent(inout) :: PLSdep(ijdim,kdim), PNSdep(ijdim,kdim)
    real(8), intent(inout) :: PLGdep(ijdim,kdim), PNGdep(ijdim,kdim)
    !+++ Melting Tendency[kg/m3/s]
    real(8), intent(in)    :: PLImlt(ijdim,kdim), PNImlt(ijdim,kdim)
    real(8), intent(in)    :: PLSmlt(ijdim,kdim), PNSmlt(ijdim,kdim)
    real(8), intent(in)    :: PLGmlt(ijdim,kdim), PNGmlt(ijdim,kdim)
    !+++
    logical, intent(in)    :: flag_history_in
    !+++ Column integrated tendency[kg/m2/s]
    real(8), intent(inout) :: sl_PLCdep(ijdim,1)
    real(8), intent(inout) :: sl_PLRdep(ijdim,1), sl_PNRdep(ijdim,1)
    !
    real(8) :: xi(ijdim,kdim)                     ! mean mass of ice particles
    real(8) :: rrhog(ijdim,kdim)                  ! 1/rhog
    real(8) :: wtem(ijdim,kdim)                   ! temperature[K]
    real(8) :: qd(ijdim,kdim)                     ! mixing ratio of dry air
    !
    real(8) :: r_cva                    ! specific heat at constant volume
    real(8) :: cpa(ijdim,kdim), r_cpa   ! specific heat at constant pressure
    real(8) :: qsw(ijdim,kdim), r_qsw   ! saturated mixing ratio for liquid
    real(8) :: qsi(ijdim,kdim), r_qsi   ! saturated mixing ratio for solid
    real(8) :: dqswdtem_rho(ijdim,kdim) ! (dqsw/dtem)_rho
    real(8) :: dqsidtem_rho(ijdim,kdim) ! (dqsi/dtem)_rho
    real(8) :: dqswdtem_pre(ijdim,kdim) ! (dqsw/dtem)_pre
    real(8) :: dqsidtem_pre(ijdim,kdim) ! (dqsi/dtem)_pre
    real(8) :: dqswdpre_tem(ijdim,kdim) ! (dqsw/dpre)_tem
    real(8) :: dqsidpre_tem(ijdim,kdim) ! (dqsi/dpre)_tem
    !
    real(8) :: w(ijdim,kdim)                     ! vetical velocity[m/s]
    real(8) :: Acnd                              ! Pdynliq + Bergeron-Findeisen
    real(8) :: Adep                              ! Pdyndep + Bergeron-Findeisen
    real(8) :: aliqliq, asolliq
    real(8) :: aliqsol, asolsol
    real(8) :: Pdynliq                           ! production term of ssw by vertical motion
    real(8) :: Pdynsol                           ! production term of ssi by vertical motion
    real(8) :: Pradliq                           ! production term of ssw by radiation
    real(8) :: Pradsol                           ! production term of ssi by radiation
    real(8) :: taucnd(ijdim,kdim),   r_taucnd    ! time scale of ssw change by MP
    real(8) :: taudep(ijdim,kdim),   r_taudep    ! time scale of ssi change by MP
    real(8) :: taucnd_c(ijdim,kdim), r_taucnd_c  ! by cloud
    real(8) :: taucnd_r(ijdim,kdim), r_taucnd_r  ! by rain 
    real(8) :: taudep_i(ijdim,kdim), r_taudep_i  ! by ice  
    real(8) :: taudep_s(ijdim,kdim), r_taudep_s  ! by snow 
    real(8) :: taudep_g(ijdim,kdim), r_taudep_g  ! by graupel
    ! alternative tendency through changing ssw and ssi
    real(8) :: PNCdep(ijdim,kdim) ! [Add] 11/08/30 T.Mitsui
    real(8) :: PLCdep_alt(ijdim,kdim)
    real(8) :: PLRdep_alt(ijdim,kdim), PNRdep_alt(ijdim,kdim)
    real(8) :: PLIdep_alt(ijdim,kdim), PNIdep_alt(ijdim,kdim)
    real(8) :: PLSdep_alt(ijdim,kdim), PNSdep_alt(ijdim,kdim)
    real(8) :: PLGdep_alt(ijdim,kdim), PNGdep_alt(ijdim,kdim)
    real(8) :: PLR2NR, PLI2NI, PLS2NS, PLG2NG
    real(8) :: coef_a_cnd, coef_b_cnd
    real(8) :: coef_a_dep, coef_b_dep
    !
    real(8) :: frz_dqc, frz_dnc
    real(8) :: frz_dqr, frz_dnr
    real(8) :: mlt_dqi, mlt_dni
    real(8) :: mlt_dqs, mlt_dns
    real(8) :: mlt_dqg, mlt_dng
    real(8) :: dep_dqi, dep_dni
    real(8) :: dep_dqs, dep_dns
    real(8) :: dep_dqg, dep_dng
    real(8) :: dep_dqr, dep_dnr
    real(8) :: dep_dqc, dep_dnc   ! 11/08/30 [Add] T.Mitsui, dep_dnc
    real(8) :: r_xc_ccn, r_xi_ccn ! 11/08/30 [Add] T.Mitsui
    !
    real(8) :: drhogqv
    real(8) :: drhogqc, drhogqr, drhogqi, drhogqs, drhogqg
    real(8) :: drhognc, drhognr, drhogni, drhogns, drhogng
    ! 
    real(8) :: fac1, fac2, fac3, fac4, fac5, fac6, fac7
    real(8) :: r_rvaptem        ! 1/(Rvap*tem)
    real(8) :: pv               ! vapor pressure 
    real(8) :: lvsw, lvsi       ! saturated vapor density
    real(8) :: dlvsw, dlvsi     ! 
    ! [Add] 11/08/30 T.Mitsui
    real(8) :: dcnd, ddep       ! total cndensation/deposition
    real(8) :: uplim_cnd        ! upper limit of condensation
    real(8) :: lowlim_cnd       ! lower limit of evaporation
    ! [Add] 11/08/30 T.Mitsui
    real(8) :: uplim_dep        ! upper limit of condensation
    real(8) :: lowlim_dep       ! lower limit of evaporation
    real(8) :: ssw, ssi         ! supersaturation ratio
    real(8) :: r_esw, r_esi     ! 1/esw, 1/esi
    real(8) :: r_lvsw, r_lvsi   ! 1/(lvsw*ssw), 1/(lvsi*ssi)
    real(8) :: evap_max         !
    real(8) :: r_dt             ! 1/dt
    real(8) :: ssw_o, ssi_o
    real(8) :: dt_dyn
    real(8) :: dt_mp
    !
    real(8)       :: tem_lh(ijdim,kdim)
    real(8)       :: dtemdt_lh(ijdim,kdim)
    real(8)       :: PLCdep_old(ijdim,kdim),PLRdep_old(ijdim,kdim)
    real(8)       :: PLIdep_old(ijdim,kdim),PLSdep_old(ijdim,kdim),PLGdep_old(ijdim,kdim)
    real(8)       :: PNRdep_old(ijdim,kdim)
    real(8)       :: PNIdep_old(ijdim,kdim),PNSdep_old(ijdim,kdim),PNGdep_old(ijdim,kdim)
    real(8)       :: PLChom_old(ijdim,kdim),PNChom_old(ijdim,kdim)
    real(8)       :: PLChet_old(ijdim,kdim),PNChet_old(ijdim,kdim)
    real(8)       :: PLRhet_old(ijdim,kdim),PNRhet_old(ijdim,kdim)
    real(8)       :: rhoge_old(ijdim,kdim), rhogq_old(ijdim,kdim,nqmax)
    real(8)       :: tem_old(ijdim,kdim), pre_old(ijdim,kdim)
    real(8), save :: fac_cndc        = 1.d0
    real(8), save :: fac_cndc_lh     = 1.d0
    logical, save :: opt_fix_taucnd_c=.false.
    logical, save :: opt_fix_lhcnd_c =.false.
    logical, save :: flag_first      =.true.
    integer :: iumax
    !
    namelist /nm_mp_ndw6_condensation/ &
         opt_fix_taucnd_c, opt_fix_lhcnd_c, fac_cndc, fac_cndc_lh

    real(8) :: fac_cndc_wrk
    !
    real(8), parameter :: tau100day   = 1.d7
    real(8), parameter :: r_tau100day = 1.d-7
    real(8), parameter :: eps=1.d-30
    !
    integer :: ij,k,iu
    !    
    ! [Add] 11/08/30 T.Mitsui
    if( flag_first )then
       flag_first = .false.
       rewind(IO_FID_CONF)
       read  (IO_FID_CONF,nml=nm_mp_ndw6_condensation, end=100)
100    write (IO_FID_LOG,nml=nm_mp_ndw6_condensation)
    end if
    !
    dt_dyn     = dt*ntmax
    dt_mp      = dt*(ntdiv-1)
    !
    r_dt       = 1.d0/dt
    rrhog(:,:) = 1.d0/rhog(:,:)
    ! Temperature lower limit is only used for saturation condition.
    ! On the other hand original "tem" is used for calculation of latent heat or energy equation.
    wtem(:,:)  = max( tem(:,:), tem_min ) 
    !
    w(:,1:kmin)    =0.D0
    w(:,kmax:kdim) =0.D0
    do k=kmin,kmax
       do ij=1,ijdim
          ! [Add] 11/08/30 T.Mitsui
          if( z(k) <= 25000.D0 )then
             w(ij,k) = 0.5d0*(wh(ij,k) + wh(ij,k+1))
          else
             w(ij,k) = 0.D0
          end if
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
    call moist_qsat_water    ( wtem, pre, qsw )
    call moist_qsat_ice      ( wtem, pre, qsi )
    call moist_dqsw_dtem_rho ( wtem, rho, dqswdtem_rho )
    call moist_dqsi_dtem_rho ( wtem, rho, dqsidtem_rho )
    call moist_dqsw_dtem_dpre( wtem, pre, dqswdtem_pre, dqswdpre_tem )
    call moist_dqsi_dtem_dpre( wtem, pre, dqsidtem_pre, dqsidpre_tem )
    ! 
    do k=1, kdim
       do ij=1, ijdim
          if( pre(ij,k) < esw(ij,k)+1.d-10 )then
             qsw(ij,k) = 1.d0
             dqswdtem_rho(ij,k) = 0.D0
             dqswdtem_pre(ij,k) = 0.D0
             dqswdpre_tem(ij,k) = 0.D0
          end if
          if( pre(ij,k) < esi(ij,k)+1.d-10 )then
             qsi(ij,k) = 1.d0
             dqsidtem_rho(ij,k) = 0.D0
             dqsidtem_pre(ij,k) = 0.D0
             dqsidpre_tem(ij,k) = 0.D0
          end if
       end do
    end do
    !
    ! taucnd, taudep
    taucnd_c(:,:)   = tau100day
    taucnd_r(:,:)   = tau100day
    taudep_i(:,:)   = tau100day
    taudep_s(:,:)   = tau100day
    taudep_g(:,:)   = tau100day
    taucnd(:,:)     = tau100day
    taudep(:,:)     = tau100day
    PLCdep_alt(:,:) = 0.D0 ! 09/08/18
    PLRdep_alt(:,:) = 0.D0
    PLIdep_alt(:,:) = 0.D0
    PLSdep_alt(:,:) = 0.D0
    PLGdep_alt(:,:) = 0.D0
    PNRdep_alt(:,:) = 0.D0
    PNIdep_alt(:,:) = 0.D0
    PNSdep_alt(:,:) = 0.D0
    PNGdep_alt(:,:) = 0.D0
    !
    if( opt_fix_lhcnd_c )then
       iumax=2
    else
       iumax=1
    end if
    !
    ! [Add] 11/08/30 T.Mitsui
    PLCdep_old(:,:) = PLCdep(:,:)
    PLRdep_old(:,:) = PLRdep(:,:)
    PLIdep_old(:,:) = PLIdep(:,:)
    PLSdep_old(:,:) = PLSdep(:,:)
    PLGdep_old(:,:) = PLGdep(:,:)
    !
    PNRdep_old(:,:) = PNRdep(:,:)
    PNIdep_old(:,:) = PNIdep(:,:)
    PNSdep_old(:,:) = PNSdep(:,:)
    PNGdep_old(:,:) = PNGdep(:,:)
    !
    PLChom_old(:,:) = PLChom(:,:)
    PLChet_old(:,:) = PLChet(:,:)
    PNChom_old(:,:) = PNChom(:,:)
    PNChet_old(:,:) = PNChet(:,:)
    PLRhet_old(:,:) = PLRhet(:,:)
    PNRhet_old(:,:) = PNRhet(:,:)
    rhoge_old(:,:)  = rhoge(:,:)
    rhogq_old(:,:,:)= rhogq(:,:,:)
    tem_old(:,:)    = tem(:,:)
    pre_old(:,:)    = pre(:,:)
    do iu=1, iumax
       if( iu==2 )then
          PLCdep(:,:) = PLCdep_old(:,:)
          PLRdep(:,:) = PLRdep_old(:,:)
          PLIdep(:,:) = PLIdep_old(:,:)
          PLSdep(:,:) = PLSdep_old(:,:)
          PLGdep(:,:) = PLGdep_old(:,:)
          !
          PNRdep(:,:) = PNRdep_old(:,:)
          PNIdep(:,:) = PNIdep_old(:,:)
          PNSdep(:,:) = PNSdep_old(:,:)
          PNGdep(:,:) = PNGdep_old(:,:)
          PLChom(:,:) = PLChom_old(:,:)
          PLChet(:,:) = PLChet_old(:,:)
          PNChom(:,:) = PNChom_old(:,:)
          PNChet(:,:) = PNChet_old(:,:)
          PLRhet(:,:) = PLRhet_old(:,:)
          PNRhet(:,:) = PNRhet_old(:,:)
          rhoge(:,:)  = rhoge_old(:,:)
          rhogq(:,:,:)= rhogq_old(:,:,:)
          tem_lh(:,:) = tem(:,:)
          tem(:,:)    = tem_old(:,:)
          pre(:,:)    = pre_old(:,:)
          !
       end if
       !
       if( opt_fix_taucnd_c )then
          fac_cndc_wrk = fac_cndc**(1.d0-b_m(I_QC))
          PLCdep(:,:)  = PLCdep_old(:,:)*fac_cndc_wrk
          write(IO_FID_LOG,*) "taucnd:fac_cndc_wrk=",fac_cndc_wrk
       end if
       !
       if( opt_fix_lhcnd_c .and. iu==1 )then
          fac_cndc_wrk = fac_cndc_lh**(1.d0-b_m(I_QC))
          PLCdep(:,:)  = PLCdep_old(:,:)*fac_cndc_wrk
          write(IO_FID_LOG,*) "lhcnd:fac_cndc_wrk=",fac_cndc_wrk          
       end if
       !
       do k=kmin, kmax
          do ij=1, ijdim
             r_rvaptem        = 1.d0/(CNST_RVAP*wtem(ij,k))
             lvsw             = esw(ij,k)*r_rvaptem        ! rho=p/(Rv*T)
             lvsi             = esi(ij,k)*r_rvaptem        ! 
             pv               = lv(ij,k)*CNST_RVAP*tem(ij,k)
             r_esw            = 1.d0/esw(ij,k)
             r_esi            = 1.d0/esi(ij,k)
             ssw              = pv*r_esw - 1.0d0 
             ssi              = pv*r_esi - 1.0d0 
             r_lvsw           = 1.d0/lvsw
             r_lvsi           = 1.d0/lvsi
             r_taucnd_c       = PLCdep(ij,k)*r_lvsw
             r_taucnd_r       = PLRdep(ij,k)*r_lvsw
             r_taudep_i       = PLIdep(ij,k)*r_lvsi
             r_taudep_s       = PLSdep(ij,k)*r_lvsi
             r_taudep_g       = PLGdep(ij,k)*r_lvsi
             taucnd_c(ij,k)   = 1.d0/(r_taucnd_c+r_tau100day)
             taucnd_r(ij,k)   = 1.d0/(r_taucnd_r+r_tau100day)
             taudep_i(ij,k)   = 1.d0/(r_taudep_i+r_tau100day)
             taudep_s(ij,k)   = 1.d0/(r_taudep_s+r_tau100day)
             taudep_g(ij,k)   = 1.d0/(r_taudep_g+r_tau100day)
             !
             r_cva            = 1.d0/cva(ij,k)
             r_cpa            = 1.d0/cpa(ij,k)
             ! Coefficient of latent heat release for ssw change by PLCdep and PLRdep
             aliqliq          = 1.d0 &
                  + r_cva*( CNST_LH00              + (CNST_CVV-CNST_CL)*tem(ij,k) )*dqswdtem_rho(ij,k)
             ! Coefficient of latent heat release for ssw change by PLIdep, PLSdep and PLGdep
             asolliq          = 1.d0 &
                  + r_cva*( CNST_LH00 + CNST_LHF00 + (CNST_CVV-CNST_CI)*tem(ij,k) )*dqswdtem_rho(ij,k)
             ! Coefficient of latent heat release for ssi change by PLCdep and PLRdep
             aliqsol          = 1.d0 &
                  + r_cva*( CNST_LH00              + (CNST_CVV-CNST_CL)*tem(ij,k) )*dqsidtem_rho(ij,k)
             ! Coefficient of latent heat release for ssi change by PLIdep, PLSdep and PLGdep
             asolsol          = 1.d0 &
                  + r_cva*( CNST_LH00 + CNST_LHF00 + (CNST_CVV-CNST_CI)*tem(ij,k) )*dqsidtem_rho(ij,k)
             Pdynliq          = w(ij,k)*CNST_EGRAV * ( r_cpa*dqswdtem_pre(ij,k) + rho(ij,k)*dqswdpre_tem(ij,k) )
             Pdynsol          = w(ij,k)*CNST_EGRAV * ( r_cpa*dqsidtem_pre(ij,k) + rho(ij,k)*dqsidpre_tem(ij,k) )  
             Pradliq          = -dTdt_rad(ij,k)    * dqswdtem_rho(ij,k)
             Pradsol          = -dTdt_rad(ij,k)    * dqsidtem_rho(ij,k)
             !
             r_qsw            = 1.d0/qsw(ij,k)
             r_qsi            = 1.d0/qsi(ij,k)
             ssw_o            = ssw - Pdynliq*r_qsw*(dt_dyn-dt_mp) + Pradliq*r_qsw*dt_mp
             ssi_o            = ssi - Pdynsol*r_qsi*(dt_dyn-dt_mp) + Pradsol*r_qsi*dt_mp
             !
             Acnd             = Pdynliq + Pradliq &
                  - ( r_taudep_i+r_taudep_s+r_taudep_g ) * ( qsw(ij,k) - qsi(ij,k) )
             Adep             = Pdynsol + Pradsol &
                  + ( r_taucnd_c+r_taucnd_r )            * ( qsw(ij,k) - qsi(ij,k) )
             r_taucnd         = &
                  + aliqliq*( r_taucnd_c+r_taucnd_r ) &
                  + asolliq*( r_taudep_i+r_taudep_s+r_taudep_g )
             r_taudep         = &
                  + aliqsol*( r_taucnd_c+r_taucnd_r )&
                  + asolsol*( r_taudep_i+r_taudep_s+r_taudep_g )
             !
!!$          uplim_cnd        = max( (lv(ij,k) - lvsw)*r_dt, 0.D0 )
!!$          lowlim_cnd       = min( (lv(ij,k) - lvsw)*r_dt, 0.D0 )
             uplim_cnd        = max( rho(ij,k)*ssw_o*qsw(ij,k)*r_dt, 0.D0 )
             lowlim_cnd       = min( rho(ij,k)*ssw_o*qsw(ij,k)*r_dt, 0.D0 )
             ! [Mod] 11/08/30 T.Mitsui
!!$             if( r_taucnd < 1.d-30 )then ! condensation is almost negligible
             if( r_taudep < r_tau100day )then
                taucnd(ij,k)     = tau100day
                PLCdep_alt(ij,k) = max(lowlim_cnd, min(uplim_cnd, PLCdep(ij,k)*ssw_o ))
                PLRdep_alt(ij,k) = max(lowlim_cnd, min(uplim_cnd, PLRdep(ij,k)*ssw_o ))
                PNRdep_alt(ij,k) = min(0.D0, PNRdep(ij,k)*ssw_o )
                PLR2NR           = 0.D0
             else
                taucnd(ij,k)     = 1.d0/r_taucnd
                ! Production term for liquid water content
                coef_a_cnd       = rho(ij,k)*Acnd*taucnd(ij,k)
                coef_b_cnd       = rho(ij,k)*taucnd(ij,k)*r_dt*(ssw_o*qsw(ij,k)-Acnd*taucnd(ij,k)) * ( exp(-dt*r_taucnd) - 1.d0 )
                PLCdep_alt(ij,k) = coef_a_cnd*r_taucnd_c - coef_b_cnd*r_taucnd_c 
                PLRdep_alt(ij,k) = coef_a_cnd*r_taucnd_r - coef_b_cnd*r_taucnd_r 
                PLR2NR           = PNRdep(ij,k)/(PLRdep(ij,k)+1.d-30)
                PNRdep_alt(ij,k) = min(0.D0, PLRdep_alt(ij,k)*PLR2NR )
             end if
             !
             uplim_dep        = max( rho(ij,k)*ssi_o*qsi(ij,k)*r_dt, 0.D0 )
             lowlim_dep       = min( rho(ij,k)*ssi_o*qsi(ij,k)*r_dt, 0.D0 )
             ! [Mod] 11/08/30 T.Mitsui
!!$             if( r_taudep < 1.d-30 )then
             if( r_taudep < r_tau100day )then
                taudep(ij,k)     = tau100day
                ! [Mod] 11/08/30 T.Mitsui, add filter for deposition by solid particle
!!$             PLIdep_alt(ij,k) = PLIdep(ij,k)*ssi_o
!!$             PLSdep_alt(ij,k) = PLSdep(ij,k)*ssi_o
!!$             PLGdep_alt(ij,k) = PLGdep(ij,k)*ssi_o
                PLIdep_alt(ij,k) = max(lowlim_dep, min(uplim_dep, PLIdep(ij,k)*ssi_o ))
                PLSdep_alt(ij,k) = max(lowlim_dep, min(uplim_dep, PLSdep(ij,k)*ssi_o ))
                PLGdep_alt(ij,k) = max(lowlim_dep, min(uplim_dep, PLGdep(ij,k)*ssi_o ))
                PNIdep_alt(ij,k) = min(0.D0, PNIdep(ij,k)*ssi_o )
                PNSdep_alt(ij,k) = min(0.D0, PNSdep(ij,k)*ssi_o )
                PNGdep_alt(ij,k) = min(0.D0, PNGdep(ij,k)*ssi_o ) 
                PLI2NI           = 0.D0 
                PLS2NS           = 0.D0 
                PLG2NG           = 0.D0 
             else
                taudep(ij,k)     = 1.d0/r_taudep 
                ! Production term for ice water content
                coef_a_dep       = rho(ij,k)*Adep*taudep(ij,k)
                coef_b_dep       = rho(ij,k)*taudep(ij,k)*r_dt*(ssi_o*qsi(ij,k)-Adep*taudep(ij,k)) * ( exp(-dt*r_taudep) - 1.d0 )
                ! [Mod] 11/08/30 T.Mitsui, add filter for deposition by solid particle
                PLIdep_alt(ij,k) =  coef_a_dep*r_taudep_i - coef_b_dep*r_taudep_i 
                PLSdep_alt(ij,k) =  coef_a_dep*r_taudep_s - coef_b_dep*r_taudep_s 
                PLGdep_alt(ij,k) =  coef_a_dep*r_taudep_g - coef_b_dep*r_taudep_g 
                ! [Mod] 11/08/30 T.Mitsui
!!$             PLI2NI           = PNIdep(ij,k)/(PLIdep(ij,k)+1.d-30)
!!$             PLS2NS           = PNSdep(ij,k)/(PLSdep(ij,k)+1.d-30)
!!$             PLG2NG           = PNGdep(ij,k)/(PLGdep(ij,k)+1.d-30)
!!$                PLI2NI           = PNIdep(ij,k)/min(PLIdep(ij,k),-1.d-30)
!!$                PLS2NS           = PNSdep(ij,k)/min(PLSdep(ij,k),-1.d-30)
!!$                PLG2NG           = PNGdep(ij,k)/min(PLGdep(ij,k),-1.d-30)
                PLI2NI           = PNIdep(ij,k)/max(PLIdep(ij,k),1.d-30)
                PLS2NS           = PNSdep(ij,k)/max(PLSdep(ij,k),1.d-30)
                PLG2NG           = PNGdep(ij,k)/max(PLGdep(ij,k),1.d-30)
                PNIdep_alt(ij,k) = min(0.D0, PLIdep_alt(ij,k)*PLI2NI )
                PNSdep_alt(ij,k) = min(0.D0, PLSdep_alt(ij,k)*PLS2NS )
                PNGdep_alt(ij,k) = min(0.D0, PLGdep_alt(ij,k)*PLG2NG )
             end if
             !
          end do
       end do
!!$       if( iu==2.and.flag_history_in )then
!!$          call history_in( 'ml_taudep', taudep(:,:) )
!!$          call history_in( 'ml_taucnd', taucnd(:,:) )
!!$          call history_in( 'ml_taucndc', taucnd_c(:,:))
!!$          call history_in( 'ml_taucndr', taucnd_r(:,:))
!!$          call history_in( 'ml_taudepi', taudep_i(:,:))
!!$          call history_in( 'ml_taudeps', taudep_s(:,:))
!!$          call history_in( 'ml_taudepg', taudep_g(:,:))
!!$       end if
       !
       PNCdep=0.D0
       PNIdep=0.D0
       !
       r_xc_ccn=1.d0/xc_ccn
       r_xi_ccn=1.d0/xi_ccn
       do k=kmin, kmax
          do ij=1, ijdim
             if( PLCdep_alt(ij,k) < -eps )then
                PNCdep(ij,k) = min(0.D0, ((lc(ij,k)+PLCdep_alt(ij,k)*dt)*r_xc_ccn - nc(ij,k))*r_dt )
             end if
             if( PLIdep(ij,k) < -eps )then
                PNIdep(ij,k) = min(0.D0, ((li(ij,k)+PLIdep_alt(ij,k)*dt)*r_xi_ccn - ni(ij,k))*r_dt )
             end if
          end do
       end do
       !
       xi(:,1:kmin)   =xi_min
       xi(:,kmax:kdim)=xi_min
       do k=kmin,kmax
          do ij=1,ijdim
             xi(ij,k) = min(xi_max, max(xi_min, li(ij,k)/(ni(ij,k)+ni_min) ))
          end do
       end do
       !
       PLCdep(:,:) = PLCdep_alt(:,:)
       PLRdep(:,:) = PLRdep_alt(:,:)
       PLIdep(:,:) = PLIdep_alt(:,:)
       PLSdep(:,:) = PLSdep_alt(:,:)
       PLGdep(:,:) = PLGdep_alt(:,:)
       PNRdep(:,:) = PNRdep_alt(:,:)
       PNIdep(:,:) = PNIdep_alt(:,:)
       PNSdep(:,:) = PNSdep_alt(:,:)
       PNGdep(:,:) = PNGdep_alt(:,:)
       !
    !OCL NORECURRENCE
       do k=kmin, kmax
          do ij=1, ijdim
             !
             !--- evaporation/condensation, deposition/sublimation
             !
             ! [Mod] 11/08/30 T.Mitsui, add filter to avoid instability
!!$          dep_dqc = max( dt*PLCdep(ij,k), -lc(ij,k) )
!!$          dep_dqr = max( dt*PLRdep(ij,k), -lr(ij,k) )
!!$          dep_dqi = max( dt*PLIdep(ij,k), -li(ij,k) )
!!$          dep_dqs = max( dt*PLSdep(ij,k), -ls(ij,k) )
!!$          dep_dqg = max( dt*PLGdep(ij,k), -lg(ij,k) )
             r_rvaptem = 1.d0/(CNST_RVAP*wtem(ij,k))
             lvsw    = esw(ij,k)*r_rvaptem
             lvsi    = esi(ij,k)*r_rvaptem
             dlvsw   = lv(ij,k)-lvsw
             dlvsi   = lv(ij,k)-lvsi  ! limiter for esi>1.d0
             dcnd    = dt*(PLCdep(ij,k)+PLRdep(ij,k))
             ddep    = dt*(PLIdep(ij,k)+PLSdep(ij,k)+PLGdep(ij,k))
             dep_dqc = 0.D0
             dep_dqr = 0.D0
             dep_dqi = 0.D0
             dep_dqs = 0.D0
             dep_dqg = 0.D0
             !    always supersaturated
             if     ( (dcnd >  eps) .and. (dlvsw > eps) )then 
                fac1    = min(dlvsw,dcnd)/dcnd
                dep_dqc =  dt*PLCdep(ij,k)*fac1
                dep_dqr =  dt*PLRdep(ij,k)*fac1
                ! always unsaturated
             else if( (dcnd < -eps) .and. (dlvsw < -eps) )then 
                fac1    = max( dlvsw,dcnd )/dcnd
                dep_dqc = max( dt*PLCdep(ij,k)*fac1, -lc(ij,k) )
                dep_dqr = max( dt*PLRdep(ij,k)*fac1, -lr(ij,k) )
             else
                ! partially unsaturated during timestep
                fac1    = 1.d0
                dep_dqc =  dt*PLCdep(ij,k)
                dep_dqr =  dt*PLRdep(ij,k)
             end if
             !
             fac2    = 1.d0
             !    always supersaturated
             if      ( (ddep >  eps) .and. (dlvsi > eps) )then 
                fac2    = min(dlvsi,ddep)/ddep
                dep_dqi = dt*PLIdep(ij,k)*fac2
                dep_dqs = dt*PLSdep(ij,k)*fac2
                dep_dqg = dt*PLGdep(ij,k)*fac2
                ! always unsaturated
             else if ( (ddep < -eps) .and. (dlvsi < -eps) )then 
                fac2    = max(dlvsi,ddep)/ddep
                dep_dqi = max(dt*PLIdep(ij,k)*fac2, -li(ij,k) )
                dep_dqs = max(dt*PLSdep(ij,k)*fac2, -ls(ij,k) )
                dep_dqg = max(dt*PLGdep(ij,k)*fac2, -lg(ij,k) )
             else
                ! partially unsaturated during timestep
                fac2    = 1.d0
                dep_dqi = dt*PLIdep(ij,k)
                dep_dqs = dt*PLSdep(ij,k)
                dep_dqg = dt*PLGdep(ij,k)
             end if
             ! evaporation always lose number(always negative).
             dep_dnc = max( dt*PNCdep(ij,k)*fac1, -nc(ij,k) ) ! ss>0 dep=0, ss<0 dep<0 ! [Add] 11/08/30 T.Mitsui
             dep_dnr = max( dt*PNRdep(ij,k)*fac1, -nr(ij,k) ) ! ss>0 dep=0, ss<0 dep<0
             dep_dni = max( dt*PNIdep(ij,k)*fac2, -ni(ij,k) ) ! ss>0 dep=0, ss<0 dep<0 
             dep_dns = max( dt*PNSdep(ij,k)*fac2, -ns(ij,k) ) ! ss>0 dep=0, ss<0 dep<0
             dep_dng = max( dt*PNGdep(ij,k)*fac2, -ng(ij,k) ) ! ss>0 dep=0, ss<0 dep<0 
             !
             !--- freezing of cloud drop
             !             
             frz_dqc = max( dt*(PLChom(ij,k)+PLChet(ij,k)), -lc(ij,k)-dep_dqc ) ! negative value
             ! [Mod] 11/08/30 T.Mitsui
!!$          frz_dnc = max( dt*(PNChom(ij,k)+PNChet(ij,k)), -nc(ij,k)         ) ! negative value 
             frz_dnc = max( dt*(PNChom(ij,k)+PNChet(ij,k)), -nc(ij,k)-dep_dnc ) ! negative value
             ! [Add] 10/08/03 T.Mitsui
             fac3    = ( frz_dqc-eps )/( dt*(PLChom(ij,k)+PLChet(ij,k))-eps )
             fac4    = ( frz_dnc-eps )/( dt*(PNChom(ij,k)+PNChet(ij,k))-eps )
             PLChom(ij,k) = fac3*PLChom(ij,k)
             PLChet(ij,k) = fac3*PLChet(ij,k)
             PNChom(ij,k) = fac4*PNChom(ij,k)
             PNChet(ij,k) = fac4*PNChet(ij,k)
             !
             !--- melting
             !
             ! ice change
             mlt_dqi = max( dt*PLImlt(ij,k), -li(ij,k)-dep_dqi )  ! negative value
             mlt_dni = max( dt*PNImlt(ij,k), -ni(ij,k)-dep_dni )  ! negative value
             ! snow change
             mlt_dqs = max( dt*PLSmlt(ij,k), -ls(ij,k)-dep_dqs )  ! negative value
             mlt_dns = max( dt*PNSmlt(ij,k), -ns(ij,k)-dep_dns )  ! negative value
             ! graupel change
             mlt_dqg = max( dt*PLGmlt(ij,k), -lg(ij,k)-dep_dqg )  ! negative value
             mlt_dng = max( dt*PNGmlt(ij,k), -ng(ij,k)-dep_dng )  ! negative value
             ! 
             !--- freezing of larger droplets
             !
             frz_dqr = max( dt*(PLRhet(ij,k)), min(0.D0, -lr(ij,k)-dep_dqr) ) ! negative value
             frz_dnr = max( dt*(PNRhet(ij,k)), min(0.D0, -nr(ij,k)-dep_dnr) ) ! negative value
             !
             ! [Add] 10/08/03 T.Mitsui
             fac5         = ( frz_dqr-eps )/( dt*PLRhet(ij,k)-eps )
             PLRhet(ij,k) = fac5*PLRhet(ij,k)
             fac6         = ( frz_dnr-eps )/( dt*PNRhet(ij,k)-eps )
             PNRhet(ij,k) = fac6*PNRhet(ij,k)
             !
             ! water vapor change
             drhogqv = -(dep_dqc+dep_dqi+dep_dqs+dep_dqg+dep_dqr)
             if ( xi(ij,k) > x_sep )then ! large ice crystals turn into rain by melting
                ! total cloud change
                drhogqc = ( frz_dqc                               + dep_dqc )
                ! [Mod] 11/08/30 T.Mitsui
!!$             drhognc = ( frz_dnc                                         )
                drhognc = ( frz_dnc                               + dep_dnc )
                ! total rain change
                drhogqr = ( frz_dqr - mlt_dqg - mlt_dqs - mlt_dqi + dep_dqr )
                drhognr = ( frz_dnr - mlt_dng - mlt_dns - mlt_dni + dep_dnr )
                !
             else
                ! total cloud change
                drhogqc = ( frz_dqc - mlt_dqi           + dep_dqc )
                ! [Mod] 11/08/30 T.Mitsui
!!$             drhognc = ( frz_dnc - mlt_dni                     )
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
             rhogq(ij,k,I_QV) = max(0.D0, rhogq(ij,k,I_QV) + drhogqv )
             !             
             rhogq(ij,k,I_QC) = max(0.D0, rhogq(ij,k,I_QC) + drhogqc )
             rhogq(ij,k,I_NC) = max(0.D0, rhogq(ij,k,I_NC) + drhognc )
             rhogq(ij,k,I_QR) = max(0.D0, rhogq(ij,k,I_QR) + drhogqr )
             rhogq(ij,k,I_NR) = max(0.D0, rhogq(ij,k,I_NR) + drhognr )
             rhogq(ij,k,I_QI) = max(0.D0, rhogq(ij,k,I_QI) + drhogqi )
             rhogq(ij,k,I_NI) = max(0.D0, rhogq(ij,k,I_NI) + drhogni )        
             rhogq(ij,k,I_QS) = max(0.D0, rhogq(ij,k,I_QS) + drhogqs )
             rhogq(ij,k,I_NS) = max(0.D0, rhogq(ij,k,I_NS) + drhogns )        
             rhogq(ij,k,I_QG) = max(0.D0, rhogq(ij,k,I_QG) + drhogqg )
             rhogq(ij,k,I_NG) = max(0.D0, rhogq(ij,k,I_NG) + drhogng )        
             !
             rhoge(ij,k) = rhoge(ij,k) &
                  - LHV * drhogqv &
                  + LHF * ( drhogqi + drhogqs + drhogqg )
             !
             !--- update mixing ratio
             q(ij,k,I_QV) = rhogq(ij,k,I_QV) * rrhog(ij,k)
             q(ij,k,I_QC) = rhogq(ij,k,I_QC) * rrhog(ij,k)
             q(ij,k,I_QR) = rhogq(ij,k,I_QR) * rrhog(ij,k)
             q(ij,k,I_QI) = rhogq(ij,k,I_QI) * rrhog(ij,k)
             q(ij,k,I_QS) = rhogq(ij,k,I_QS) * rrhog(ij,k)
             q(ij,k,I_QG) = rhogq(ij,k,I_QG) * rrhog(ij,k)
             !
             q(ij,k,I_NC) = rhogq(ij,k,I_NC) * rrhog(ij,k)
             q(ij,k,I_NR) = rhogq(ij,k,I_NR) * rrhog(ij,k)
             q(ij,k,I_NI) = rhogq(ij,k,I_NI) * rrhog(ij,k)
             q(ij,k,I_NS) = rhogq(ij,k,I_NS) * rrhog(ij,k)
             q(ij,k,I_NG) = rhogq(ij,k,I_NG) * rrhog(ij,k)
             ! 
             sl_PLCdep(ij,1) = sl_PLCdep(ij,1) + dep_dqc*Dz(k)*gsgam2(ij,k)
             sl_PLRdep(ij,1) = sl_PLRdep(ij,1) + dep_dqr*Dz(k)*gsgam2(ij,k)
             sl_PNRdep(ij,1) = sl_PNRdep(ij,1) + dep_dnr*Dz(k)*gsgam2(ij,k)
          end do
       end do
       call thrmdyn_qd(   &
            qd,           & !--- out
            q )             !--- in
       call thrmdyn_cv(   &
            cva,          & !--- out
            q,            & !--- in
            qd )            !--- in
       tem(:,:) = rhoge(:,:) / ( rhog(:,:) * cva(:,:) )
       !    
       ! [Add] 11/08/30 T.Mitsui
       if( opt_fix_lhcnd_c .and. iu==2 )then
          rhoge(:,:) = tem_lh(:,:)*rhog(:,:)*cva(:,:)
          if( flag_history_in )then
             dtemdt_lh(:,:) = (tem(:,:)-tem_lh(:,:))*r_dt
!!$             call history_in( 'ml_dTdt_lh', dtemdt_lh(:,:) )
          end if
          tem(:,:)   = tem_lh(:,:)
       end if
       pre(:,:) = rho(:,:)*( qd(:,:)*CNST_RAIR+q(:,:,I_QV)*CNST_RVAP )*tem(:,:)
    end do
    !
    return
  end subroutine update_by_phase_change
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
   subroutine negative_filter
    use mod_grid, only: &
       KA => GRID_KA, &
       IA => GRID_IA, &
       JA => GRID_JA, &
       KS => GRID_KS, &
       KE => GRID_KE, &
       IS => GRID_IS, &
       IE => GRID_IE, &
       JS => GRID_JS, &
       JE => GRID_JE
     use mod_atmos_vars, only: &
       var => atmos_var, &
       QA  => A_QA, &
       I_DENS,      &
       I_QV,        &
       I_QC,        &
       I_QR,        &
       I_QI,        &
       I_QS,        &
       I_QG,        &
       I_NC,        &
       I_NR,        &
       I_NI,        &
       I_NS,        &
       I_NG
    implicit none

    real(8) :: diffq(KA,IA,JA)
    real(8) :: r_xmin

    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    ! total hydrometeor (before correction)
    !OCL NORECURRENCE
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       diffq(k,i,j) = var(k,i,j,5+I_QV) &
                    + var(k,i,j,5+I_QC) &
                    + var(k,i,j,5+I_QR) &
                    + var(k,i,j,5+I_QI) &
                    + var(k,i,j,5+I_QS) &
                    + var(k,i,j,5+I_QG)
    enddo
    enddo
    enddo

    ! remove negative value of hydrometeor (mass, number)
    !OCL SIMD
    do iq = 1, QA
    do j  = 1, JA
    do i  = 1, IA
    do k  = 1, KA
       if ( var(k,i,j,5+iq) < 0.D0 ) then
          var(k,i,j,5+iq) = 0.D0
       endif
    enddo
    enddo
    enddo
    enddo

    ! apply correction of hydrometeor to total density
    ! [note] mass conservation is broken here to fill rounding error.
    !OCL NORECURRENCE
    do j  = 1, JA
    do i  = 1, IA
    do k  = 1, KA
       var(k,i,j,I_DENS) = var(k,i,j,I_DENS)   &
                         * ( 1.D0              &
                           + var(k,i,j,5+I_QV) &
                           + var(k,i,j,5+I_QC) &
                           + var(k,i,j,5+I_QR) &
                           + var(k,i,j,5+I_QI) &
                           + var(k,i,j,5+I_QS) &
                           + var(k,i,j,5+I_QG) &
                           - diffq(k,i,j)      ) ! after-before
    enddo
    enddo
    enddo

    ! avoid unrealistical value of number concentration 
    ! due to numerical diffusion in advection
    r_xmin = 1.D0 / xmin_filter    

    !OCL SIMD
    do j  = 1, JA
    do i  = 1, IA
    do k  = 1, KA
       if ( var(k,i,j,5+I_NC) > var(k,i,j,5+I_QC)*r_xmin ) then
          var(k,i,j,5+I_NC) = var(k,i,j,5+I_QC)*r_xmin
       endif
    enddo
    enddo
    enddo
    !OCL SIMD
    do j  = 1, JA
    do i  = 1, IA
    do k  = 1, KA
       if ( var(k,i,j,5+I_NR) > var(k,i,j,5+I_QR)*r_xmin ) then
          var(k,i,j,5+I_NR) = var(k,i,j,5+I_QR)*r_xmin
       endif
    enddo
    enddo
    enddo
    !OCL SIMD
    do j  = 1, JA
    do i  = 1, IA
    do k  = 1, KA
       if ( var(k,i,j,5+I_NI) > var(k,i,j,5+I_QI)*r_xmin ) then
          var(k,i,j,5+I_NI) = var(k,i,j,5+I_QI)*r_xmin
       endif
    enddo
    enddo
    enddo
    !OCL SIMD
    do j  = 1, JA
    do i  = 1, IA
    do k  = 1, KA
       if ( var(k,i,j,5+I_NS) > var(k,i,j,5+I_QS)*r_xmin ) then
          var(k,i,j,5+I_NS) = var(k,i,j,5+I_QS)*r_xmin
       endif
    enddo
    enddo
    enddo
    !OCL SIMD
    do j  = 1, JA
    do i  = 1, IA
    do k  = 1, KA
       if ( var(k,i,j,5+I_NG) > var(k,i,j,5+I_QG)*r_xmin ) then
          var(k,i,j,5+I_NG) = var(k,i,j,5+I_QG)*r_xmin
       endif
    enddo
    enddo
    enddo

    return
  end subroutine negative_filter


end module mod_atmos_phy_mp
!-------------------------------------------------------------------------------
