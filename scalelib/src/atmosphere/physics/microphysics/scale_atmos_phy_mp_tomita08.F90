!-------------------------------------------------------------------------------
!> module atmosphere / physics / microphysics / Tomita08
!!
!! @par Description
!!          Cloud Microphysics by Lin-type parametarization
!!          Reference: Tomita(2008)
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_phy_mp_tomita08
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_MP_tomita08_setup
  public :: ATMOS_PHY_MP_tomita08_adjustment
  public :: ATMOS_PHY_MP_tomita08_terminal_velocity
  public :: ATMOS_PHY_MP_tomita08_effective_radius
  public :: ATMOS_PHY_MP_tomita08_cloud_fraction
  public :: ATMOS_PHY_MP_tomita08_qtrc2qhyd
  public :: ATMOS_PHY_MP_tomita08_qhyd2qtrc

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, private, parameter :: QA_MP  = 6

  integer,                parameter, public :: ATMOS_PHY_MP_tomita08_ntracers = QA_MP
  integer,                parameter, public :: ATMOS_PHY_MP_tomita08_nwaters = 2
  integer,                parameter, public :: ATMOS_PHY_MP_tomita08_nices = 3
  character(len=H_SHORT), parameter, public :: ATMOS_PHY_MP_tomita08_tracer_names(QA_MP) = (/ &
       'QV', &
       'QC', &
       'QR', &
       'QI', &
       'QS', &
       'QG'  /)
  character(len=H_MID)  , parameter, public :: ATMOS_PHY_MP_tomita08_tracer_descriptions(QA_MP) = (/ &
       'Ratio of Water Vapor mass to total mass (Specific humidity)', &
       'Ratio of Cloud Water mass to total mass                    ', &
       'Ratio of Rain Water mass to total mass                     ', &
       'Ratio of Cloud Ice mass ratio to total mass                ', &
       'Ratio of Snow miass ratio to total mass                    ', &
       'Ratio of Graupel mass ratio to total mass                  '/)
  character(len=H_SHORT), parameter, public :: ATMOS_PHY_MP_tomita08_tracer_units(QA_MP) = (/ &
       'kg/kg',  &
       'kg/kg',  &
       'kg/kg',  &
       'kg/kg',  &
       'kg/kg',  &
       'kg/kg'   /)

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: MP_tomita08
  private :: MP_tomita08_BergeronParam

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,  private, parameter   :: I_QV = 1
  integer,  private, parameter   :: I_QC = 2
  integer,  private, parameter   :: I_QR = 3
  integer,  private, parameter   :: I_QI = 4
  integer,  private, parameter   :: I_QS = 5
  integer,  private, parameter   :: I_QG = 6

  integer,  private, parameter   :: I_hyd_QC =  1
  integer,  private, parameter   :: I_hyd_QR =  2
  integer,  private, parameter   :: I_hyd_QI =  3
  integer,  private, parameter   :: I_hyd_QS =  4
  integer,  private, parameter   :: I_hyd_QG =  5

  logical,  private              :: do_couple_aerosol   ! apply CCN effect?
  logical,  private              :: do_explicit_icegen  ! apply explicit ice generation?

  logical,  private              :: fixed_re  = .false. ! use ice's effective radius for snow and graupel, and set rain transparent?
  logical,  private              :: const_rec = .true.  ! use constant  effective radius for cloud water?
  logical,  private              :: nofall_qr = .false. ! surpress sedimentation of rain?
  logical,  private              :: nofall_qi = .false. ! surpress sedimentation of ice?
  logical,  private              :: nofall_qs = .false. ! surpress sedimentation of snow?
  logical,  private              :: nofall_qg = .false. ! surpress sedimentation of graupel?

  real(RP), private, parameter   :: dens00      = 1.28_RP !< standard density [kg/m3]

  ! Parameter for Marshall-Palmer distribution
  real(RP), private              :: N0r_def     = 8.E+6_RP !< intercept parameter for rain    [1/m4]
  real(RP), private              :: N0s_def     = 3.E+6_RP !< intercept parameter for snow    [1/m4]
  real(RP), private              :: N0g_def     = 4.E+6_RP !< intercept parameter for graupel [1/m4]

  real(RP), private              :: dens_s      = 100.0_RP !< density of snow    [kg/m3]
  real(RP), private              :: dens_g      = 400.0_RP !< density of graupel [kg/m3]
                                                           !   graupel : 400
                                                           !   hail    : 917
  real(RP), private              :: drag_g      = 0.6_RP   !< drag coefficient for graupel
  real(RP), private              :: re_qc       =  8.E-6_RP ! effective radius for cloud water
  real(RP), private              :: re_qi       = 40.E-6_RP ! effective radius for cloud ice
  real(RP), private              :: Cr          = 130.0_RP
  real(RP), private              :: Cs          = 4.84_RP

  ! Empirical parameter
  real(RP), private              :: Ar, As, Ag
  real(RP), private              :: Br, Bs, Bg
  real(RP), private              :: Cg
  real(RP), private              :: Dr, Ds, Dg

  ! GAMMA function
  real(RP), private              :: GAM, GAM_2, GAM_3

  real(RP), private              :: GAM_1br, GAM_2br, GAM_3br
  real(RP), private              :: GAM_3dr
  real(RP), private              :: GAM_6dr
  real(RP), private              :: GAM_1brdr
  real(RP), private              :: GAM_5dr_h

  real(RP), private              :: GAM_1bs, GAM_2bs, GAM_3bs
  real(RP), private              :: GAM_3ds
  real(RP), private              :: GAM_1bsds
  real(RP), private              :: GAM_5ds_h

  real(RP), private              :: GAM_1bg, GAM_3dg
  real(RP), private              :: GAM_1bgdg
  real(RP), private              :: GAM_5dg_h

  !---< Khairoutdinov and Kogan (2000) >---
  logical,  private              :: enable_KK2000   = .false. !< use scheme by Khairoutdinov and Kogan (2000)


  !---< Roh and Satoh (2014) >---
  logical,  private              :: enable_RS2014   = .false. !< use scheme by Roh and Satoh (2014)
  real(RP), private              :: ln10                 !< log(10)
  real(RP), private, parameter   :: coef_a01 =  5.065339_RP
  real(RP), private, parameter   :: coef_a02 = -0.062659_RP
  real(RP), private, parameter   :: coef_a03 = -3.032362_RP
  real(RP), private, parameter   :: coef_a04 =  0.029469_RP
  real(RP), private, parameter   :: coef_a05 = -0.000285_RP
  real(RP), private, parameter   :: coef_a06 =  0.31255_RP
  real(RP), private, parameter   :: coef_a07 =  0.000204_RP
  real(RP), private, parameter   :: coef_a08 =  0.003199_RP
  real(RP), private, parameter   :: coef_a09 =  0.0_RP
  real(RP), private, parameter   :: coef_a10 = -0.015952_RP

  real(RP), private, parameter   :: coef_b01 =  0.476221_RP
  real(RP), private, parameter   :: coef_b02 = -0.015896_RP
  real(RP), private, parameter   :: coef_b03 =  0.165977_RP
  real(RP), private, parameter   :: coef_b04 =  0.007468_RP
  real(RP), private, parameter   :: coef_b05 = -0.000141_RP
  real(RP), private, parameter   :: coef_b06 =  0.060366_RP
  real(RP), private, parameter   :: coef_b07 =  0.000079_RP
  real(RP), private, parameter   :: coef_b08 =  0.000594_RP
  real(RP), private, parameter   :: coef_b09 =  0.0_RP
  real(RP), private, parameter   :: coef_b10 = -0.003577_RP

  !---< Wainwright et al. (2014) >---
  logical,  private              :: enable_WDXZ2014 = .false. !< use scheme by Wainwright et al. (2014)

  ! Accretion parameter
  real(RP), private              :: Eiw         = 1.0_RP      !< collection efficiency of cloud ice for cloud water
  real(RP), private              :: Erw         = 1.0_RP      !< collection efficiency of rain    for cloud water
  real(RP), private              :: Esw         = 1.0_RP      !< collection efficiency of snow    for cloud water
  real(RP), private              :: Egw         = 1.0_RP      !< collection efficiency of graupel for cloud water
  real(RP), private              :: Eri         = 1.0_RP      !< collection efficiency of rain    for cloud ice
  real(RP), private              :: Esi         = 1.0_RP      !< collection efficiency of snow    for cloud ice
  real(RP), private              :: Egi         = 0.1_RP      !< collection efficiency of graupel for cloud ice
  real(RP), private              :: Esr         = 1.0_RP      !< collection efficiency of snow    for rain
  real(RP), private              :: Egr         = 1.0_RP      !< collection efficiency of graupel for rain
  real(RP), private              :: Egs         = 1.0_RP      !< collection efficiency of graupel for snow
  real(RP), private              :: gamma_sacr  = 25.E-3_RP   !< effect of low temperature for Esi
  real(RP), private              :: gamma_gacs  = 90.E-3_RP   !< effect of low temperature for Egs
  real(RP), private              :: mi          = 4.19E-13_RP !< mass of one cloud ice crystal [kg]

  ! Auto-conversion parameter
  real(RP), private, parameter   :: Nc_lnd      = 2000.0_RP   !< number concentration of cloud water (land)  [1/cc]
  real(RP), private, parameter   :: Nc_ocn      =   50.0_RP   !< number concentration of cloud water (ocean) [1/cc]
  real(RP), private, allocatable :: Nc_def(:,:)               !< number concentration of cloud water         [1/cc]

  real(RP), private              :: beta_saut   =  6.E-3_RP   !< auto-conversion factor beta  for ice
  real(RP), private              :: gamma_saut  = 60.E-3_RP   !< auto-conversion factor gamma for ice
  real(RP), private              :: beta_gaut   =  0.0_RP     !< auto-conversion factor beta  for snow
  real(RP), private              :: gamma_gaut  = 90.E-3_RP   !< auto-conversion factor gamma for snow
  real(RP), private              :: qicrt_saut  =  0.0_RP     !< mass ratio threshold for Psaut [kg/kg]
  real(RP), private              :: qscrt_gaut  =  6.E-4_RP   !< mass ratio threshold for Pgaut [kg/kg]

  ! Evaporation, Sublimation parameter
  real(RP), private, parameter   :: Da0         = 2.428E-2_RP !< thermal diffusion coefficient of air at 0C,1atm [J/m/s/K]
  real(RP), private, parameter   :: dDa_dT      =  7.47E-5_RP !< Coefficient of Da depending on temperature      [J/m/s/K/K]
  real(RP), private, parameter   :: Dw0         = 2.222E-5_RP !< diffusion coefficient of water vapor in the air at 0C,1atm [m2/s]
  real(RP), private, parameter   :: dDw_dT      =  1.37E-7_RP !< Coefficient of Dw depending on temperature                 [m2/s/K]
  real(RP), private, parameter   :: mu0         = 1.718E-5_RP !< kinematic viscosity of air at 0C,1atm      [m2/s*kg/m3]
  real(RP), private, parameter   :: dmu_dT      =  5.28E-8_RP !< Coefficient of mu depending on temperature [m2/s/K*kg/m3]

  real(RP), private, parameter   :: f1r         = 0.78_RP     !< ventilation factor 1 for rain
  real(RP), private, parameter   :: f2r         = 0.27_RP     !< ventilation factor 2 for rain    ( x Schumidt number? )
  real(RP), private, parameter   :: f1s         = 0.65_RP     !< ventilation factor 1 for snow
  real(RP), private, parameter   :: f2s         = 0.39_RP     !< ventilation factor 2 for snow    ( x Schumidt number? )
  real(RP), private, parameter   :: f1g         = 0.78_RP     !< ventilation factor 1 for graupel
  real(RP), private, parameter   :: f2g         = 0.27_RP     !< ventilation factor 2 for graupel ( x Schumidt number? )

  ! Freezing parameter
  real(RP), private, parameter   :: A_frz       =  0.66_RP    !< freezing factor [/K]
  real(RP), private, parameter   :: B_frz       = 100.0_RP    !< freezing factor [/m3/s]

  ! Bergeron process parameter
  real(RP), private, parameter   :: mi40        = 2.46E-10_RP !< mass              of a 40 micron ice crystal [kg]
  real(RP), private, parameter   :: mi50        = 4.80E-10_RP !< mass              of a 50 micron ice crystal [kg]
  real(RP), private, parameter   :: vti50       = 1.0_RP      !< terminal velocity of a 50 micron ice crystal [m/s]
  real(RP), private, parameter   :: Ri50        = 5.E-5_RP    !< radius            of a 50 micron ice crystal [m]

  ! Explicit ice generation
  logical,  private              :: only_liquid = .false.
  real(RP), private              :: sw_expice   = 0.0_RP
  real(RP), private, parameter   :: Nc_ihtr     = 300.0_RP  !< cloud number concentration for heterogeneous ice nucleation [1/cc]
  real(RP), private, parameter   :: Di_max      = 500.E-6_RP
  real(RP), private, parameter   :: Di_a        = 11.9_RP

  integer,  private, parameter   :: w_nmax = 49
  integer,  private, parameter   :: I_dqv_dt  =  1 !
  integer,  private, parameter   :: I_dqc_dt  =  2 !
  integer,  private, parameter   :: I_dqr_dt  =  3 !
  integer,  private, parameter   :: I_dqi_dt  =  4 !
  integer,  private, parameter   :: I_dqs_dt  =  5 !
  integer,  private, parameter   :: I_dqg_dt  =  6 !
  integer,  private, parameter   :: I_delta1  =  7 ! separation switch for r->s,g
  integer,  private, parameter   :: I_delta2  =  8 ! separation switch for s->g
  integer,  private, parameter   :: I_spsati  =  9 ! separation switch for ice sublimation
  integer,  private, parameter   :: I_iceflg  = 10 ! separation switch for T > 0
  integer,  private, parameter   :: I_RLMDr   = 11
  integer,  private, parameter   :: I_RLMDs   = 12
  integer,  private, parameter   :: I_RLMDg   = 13
  integer,  private, parameter   :: I_Piacr   = 14 ! r->s,g
  integer,  private, parameter   :: I_Psacr   = 15 ! r->s,g
  integer,  private, parameter   :: I_Praci   = 16 ! i->s,g
  integer,  private, parameter   :: I_Pigen   = 17 ! v->i
  integer,  private, parameter   :: I_Pidep   = 18 ! v->i
  integer,  private, parameter   :: I_Psdep   = 19 ! v->s
  integer,  private, parameter   :: I_Pgdep   = 20 ! v->g
  integer,  private, parameter   :: I_Praut   = 21 ! c->r
  integer,  private, parameter   :: I_Pracw   = 22 ! c->r
  integer,  private, parameter   :: I_Pihom   = 23 ! c->i
  integer,  private, parameter   :: I_Pihtr   = 24 ! c->i
  integer,  private, parameter   :: I_Psacw   = 25 ! c->s
  integer,  private, parameter   :: I_Psfw    = 26 ! c->s
  integer,  private, parameter   :: I_Pgacw   = 27 ! c->g
  integer,  private, parameter   :: I_Prevp   = 28 ! r->v
  integer,  private, parameter   :: I_Piacr_s = 29 ! r->s
  integer,  private, parameter   :: I_Psacr_s = 30 ! r->s
  integer,  private, parameter   :: I_Piacr_g = 31 ! r->g
  integer,  private, parameter   :: I_Psacr_g = 32 ! r->g
  integer,  private, parameter   :: I_Pgacr   = 33 ! r->g
  integer,  private, parameter   :: I_Pgfrz   = 34 ! r->g
  integer,  private, parameter   :: I_Pisub   = 35 ! i->v
  integer,  private, parameter   :: I_Pimlt   = 36 ! i->c
  integer,  private, parameter   :: I_Psaut   = 37 ! i->s
  integer,  private, parameter   :: I_Praci_s = 38 ! i->s
  integer,  private, parameter   :: I_Psaci   = 39 ! i->s
  integer,  private, parameter   :: I_Psfi    = 40 ! i->s
  integer,  private, parameter   :: I_Praci_g = 41 ! i->g
  integer,  private, parameter   :: I_Pgaci   = 42 ! i->g
  integer,  private, parameter   :: I_Pssub   = 43 ! s->v
  integer,  private, parameter   :: I_Psmlt   = 44 ! s->r
  integer,  private, parameter   :: I_Pgaut   = 45 ! s->g
  integer,  private, parameter   :: I_Pracs   = 46 ! s->g
  integer,  private, parameter   :: I_Pgacs   = 47 ! s->g
  integer,  private, parameter   :: I_Pgsub   = 48 ! g->v
  integer,  private, parameter   :: I_Pgmlt   = 49 ! g->r

  character(len=H_SHORT), private :: w_name(w_nmax)

  data w_name / 'dqv_dt ', &
                'dqc_dt ', &
                'dqr_dt ', &
                'dqi_dt ', &
                'dqs_dt ', &
                'dqg_dt ', &
                'delta1 ', &
                'delta2 ', &
                'spsati ', &
                'iceflg ', &
                'RLMDr  ', &
                'RLMDs  ', &
                'RLMDg  ', &
                'Piacr  ', &
                'Psacr  ', &
                'Praci  ', &
                'Pigen  ', &
                'Pidep  ', &
                'Psdep  ', &
                'Pgdep  ', &
                'Praut  ', &
                'Pracw  ', &
                'Pihom  ', &
                'Pihtr  ', &
                'Psacw  ', &
                'Psfw   ', &
                'Pgacw  ', &
                'Prevp  ', &
                'Piacr_s', &
                'Psacr_s', &
                'Piacr_g', &
                'Psacr_g', &
                'Pgacr  ', &
                'Pgfrz  ', &
                'Pisub  ', &
                'Pimlt  ', &
                'Psaut  ', &
                'Praci_s', &
                'Psaci  ', &
                'Psfi   ', &
                'Praci_g', &
                'Pgaci  ', &
                'Pssub  ', &
                'Psmlt  ', &
                'Pgaut  ', &
                'Pracs  ', &
                'Pgacs  ', &
                'Pgsub  ', &
                'Pgmlt  '  /

  real(RP), private, allocatable :: w3d(:,:,:,:) !< for history output
  integer,  private              :: HIST_id(w_nmax)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> ATMOS_PHY_MP_tomita08_setup
  !! Setup
  !<
  subroutine ATMOS_PHY_MP_tomita08_setup( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE )
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       PI     => CONST_PI,    &
       GRAV   => CONST_GRAV,  &
       dens_w => CONST_DWATR, &
       dens_i => CONST_DICE
    use scale_specfunc, only: &
       SF_gamma
    use scale_file_history, only: &
       FILE_HISTORY_reg
    implicit none

    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP) :: autoconv_nc     = Nc_ocn  !< number concentration of cloud water [1/cc]

    namelist / PARAM_ATMOS_PHY_MP_TOMITA08 / &
       do_couple_aerosol, &
       do_explicit_icegen, &
       autoconv_nc,     &
       enable_KK2000,   &
       enable_RS2014,   &
       enable_WDXZ2014, &
       N0r_def,         &
       N0s_def,         &
       N0g_def,         &
       dens_s,          &
       dens_g,          &
       re_qc,           &
       re_qi,           &
       drag_g,          &
       Cr,              &
       Cs,              &
       Eiw,             &
       Erw,             &
       Esw,             &
       Egw,             &
       Eri,             &
       Esi,             &
       Egi,             &
       Esr,             &
       Egr,             &
       Egs,             &
       gamma_sacr,      &
       gamma_gacs,      &
       mi,              &
       beta_saut,       &
       gamma_saut,      &
       qicrt_saut,      &
       beta_gaut,       &
       gamma_gaut,      &
       qscrt_gaut,      &
       fixed_re,        &
       const_rec,       &
       nofall_qr,       &
       nofall_qi,       &
       nofall_qs,       &
       nofall_qg

    real(RP), parameter :: max_term_vel = 10.0_RP  !-- terminal velocity for calculate dt of sedimentation

    integer  :: ierr
    integer  :: i, j, ip
    !---------------------------------------------------------------------------


    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_MP_tomita08_setup",*) 'Setup'
    LOG_INFO("ATMOS_PHY_MP_tomita08_setup",*) 'Tomita (2008) 1-moment bulk 6 category'

    allocate( w3d(KA,IA,JA,w_nmax) )
    w3d(:,:,:,:) = 0.0_RP

    allocate( Nc_def(IA,JA) )


    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_MP_TOMITA08,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_PHY_MP_tomita08_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_MP_tomita08_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_MP_TOMITA08. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_MP_TOMITA08)

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_MP_tomita08_setup",*) 'density of the snow    [kg/m3] : ', dens_s
    LOG_INFO("ATMOS_PHY_MP_tomita08_setup",*) 'density of the graupel [kg/m3] : ', dens_g
    LOG_INFO("ATMOS_PHY_MP_tomita08_setup",*) 'Nc for auto-conversion [num/m3]: ', autoconv_nc
    LOG_INFO("ATMOS_PHY_MP_tomita08_setup",*) 'Use k-k  scheme?               : ', enable_KK2000
    LOG_INFO("ATMOS_PHY_MP_tomita08_setup",*) 'Use Roh  scheme?               : ', enable_RS2014
    LOG_INFO("ATMOS_PHY_MP_tomita08_setup",*) 'Use WDXZ scheme?               : ', enable_WDXZ2014
    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_MP_tomita08_setup",*) 'Use effective radius of ice for snow and graupel,'
    LOG_INFO("ATMOS_PHY_MP_tomita08_setup",*) '    and set rain transparent?          : ', fixed_re
    LOG_INFO("ATMOS_PHY_MP_tomita08_setup",*) 'Density of the ice is used for the calculation of '
    LOG_INFO("ATMOS_PHY_MP_tomita08_setup",*) '    optically effective volume of snow and graupel.'
    LOG_INFO("ATMOS_PHY_MP_tomita08_setup",*) 'Surpress sedimentation of rain?    : ', nofall_qr
    LOG_INFO("ATMOS_PHY_MP_tomita08_setup",*) 'Surpress sedimentation of ice?     : ', nofall_qi
    LOG_INFO("ATMOS_PHY_MP_tomita08_setup",*) 'Surpress sedimentation of snow?    : ', nofall_qs
    LOG_INFO("ATMOS_PHY_MP_tomita08_setup",*) 'Surpress sedimentation of graupel? : ', nofall_qg
    LOG_INFO("ATMOS_PHY_MP_tomita08_setup",*) 'Enable explicit ice generation?    : ', do_explicit_icegen
    LOG_NEWLINE

    do j = JS, JE
    do i = IS, IE
       Nc_def(i,j) = autoconv_nc
    end do
    end do

    !--- empirical coefficients A, B, C, D
    Ar = PI * dens_w / 6.0_RP
    As = PI * dens_s / 6.0_RP
    Ag = PI * dens_g / 6.0_RP

    Br = 3.0_RP
    Bs = 3.0_RP
    Bg = 3.0_RP

    Cg = sqrt( ( 4.0_RP * dens_g * GRAV ) / ( 3.0_RP * dens00 * drag_g ) )

    Dr = 0.50_RP
    Ds = 0.25_RP
    Dg = 0.50_RP

    if ( enable_RS2014 ) then ! overwrite parameters
       do_explicit_icegen = .true.

       N0g_def   = 4.E+8_RP
       As        = 0.069_RP
       Bs        = 2.0_RP
       Esi       = 0.25_RP
       Egi       = 0.0_RP
       Egs       = 0.0_RP
    endif

    if ( do_explicit_icegen ) then
       only_liquid = .true.
       sw_expice   = 1.0_RP
    else
       only_liquid = .false.
       sw_expice   = 0.0_RP
    endif

    GAM       = 1.0_RP ! =0!
    GAM_2     = 1.0_RP ! =1!
    GAM_3     = 2.0_RP ! =2!

    GAM_1br   = SF_gamma( 1.0_RP + Br ) ! = 4!
    GAM_2br   = SF_gamma( 2.0_RP + Br ) ! = 5!
    GAM_3br   = SF_gamma( 3.0_RP + Br ) ! = 6!
    GAM_3dr   = SF_gamma( 3.0_RP + Dr )
    GAM_6dr   = SF_gamma( 6.0_RP + Dr )
    GAM_1brdr = SF_gamma( 1.0_RP + Br + Dr )
    GAM_5dr_h = SF_gamma( 0.5_RP * (5.0_RP+Dr) )

    GAM_1bs   = SF_gamma( 1.0_RP + Bs ) ! = 4!
    GAM_2bs   = SF_gamma( 2.0_RP + Bs ) ! = 5!
    GAM_3bs   = SF_gamma( 3.0_RP + Bs ) ! = 6!
    GAM_3ds   = SF_gamma( 3.0_RP + Ds )
    GAM_1bsds = SF_gamma( 1.0_RP + Bs + Ds )
    GAM_5ds_h = SF_gamma( 0.5_RP * (5.0_RP+Ds) )

    GAM_1bg   = SF_gamma( 1.0_RP + Bg ) ! = 4!
    GAM_3dg   = SF_gamma( 3.0_RP + Dg )
    GAM_1bgdg = SF_gamma( 1.0_RP + Bg + Dg)
    GAM_5dg_h = SF_gamma( 0.5_RP * (5.0_RP+Dg) )

    ln10 = log(10.0_RP)

    ! history
    do ip = 1, w_nmax
       call FILE_HISTORY_reg( w_name(ip), 'individual tendency term in tomita08', 'kg/kg/s', & ! [IN]
                              hist_id(ip)                                                    ) ! [OUT]
    end do

    return
  end subroutine ATMOS_PHY_MP_tomita08_setup

  !-----------------------------------------------------------------------------
  !> ATMOS_PHY_MP_tomita08_adjustment
  !! calculate state after saturation process
  !<
  subroutine ATMOS_PHY_MP_tomita08_adjustment( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       DENS, PRES,               &
       CCN,                      &
       dt,                       &
       TEMP, QTRC, CPtot, CVtot, &
       RHOE_t, EVAPORATE         )
    use scale_const, only: &
       DWATR => CONST_DWATR, &
       PI    => CONST_PI
    use scale_file_history, only: &
       FILE_HISTORY_in
    use scale_atmos_phy_mp_common, only: &
       MP_saturation_adjustment => ATMOS_PHY_MP_saturation_adjustment
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in) :: DENS (KA,IA,JA)
    real(RP), intent(in) :: PRES (KA,IA,JA)
    real(RP), intent(in) :: CCN  (KA,IA,JA)
    real(DP), intent(in) :: dt

    real(RP), intent(inout) :: TEMP(KA,IA,JA)
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA_MP)
    real(RP), intent(inout) :: CPtot(KA,IA,JA)
    real(RP), intent(inout) :: CVtot(KA,IA,JA)

    real(RP), intent(out) :: RHOE_t   (KA,IA,JA)
    real(RP), intent(out) :: EVAPORATE(KA,IA,JA)   ! number of evaporated cloud [/m3/s]

    real(RP) :: RHOE_d_sat(KA,IA,JA)

    real(RP) :: QC_t_sat(KA,IA,JA)
    real(RP) :: QI_t_sat(KA,IA,JA)

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    LOG_PROGRESS(*) 'atmosphere / physics / microphysics / Tomita08'

    !##### MP Main #####
    call MP_tomita08( &
         KA, KS, KE, IA, IS, IE, JA, JS, JE, &
         DENS(:,:,:), PRES(:,:,:), CCN(:,:,:), & ! [IN]
         dt,                                   & ! [IN]
         TEMP(:,:,:), QTRC(:,:,:,:),           & ! [INOUT]
         CPtot(:,:,:), CVtot(:,:,:),           & ! [INOUT]
         RHOE_t(:,:,:)                         ) ! [OUT]

    ! save value before saturation adjustment
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       QC_t_sat(k,i,j) = QTRC(k,i,j,I_QC)
       QI_t_sat(k,i,j) = QTRC(k,i,j,I_QI)
    enddo
    enddo
    enddo

    call MP_saturation_adjustment( &
         KA, KS, KE, IA, IS, IE, JA, JS, JE, &
         DENS(:,:,:),                        & ! [IN]
         only_liquid,                        & ! [IN]
         TEMP(:,:,:),                        & ! [INOUT]
         QTRC(:,:,:,I_QV),                   & ! [INOUT]
         QTRC(:,:,:,I_QC), QTRC(:,:,:,I_QI), & ! [INOUT]
         CPtot(:,:,:), CVtot(:,:,:),         & ! [INOUT]
         RHOE_d_sat(:,:,:)                   ) ! [OUT]

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       RHOE_t(k,i,j) = RHOE_t(k,i,j) + RHOE_d_sat(k,i,j) / dt
    enddo
    enddo
    enddo
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       QC_t_sat(k,i,j) = ( QTRC(k,i,j,I_QC) - QC_t_sat(k,i,j) ) / dt
    enddo
    enddo
    enddo
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       QI_t_sat(k,i,j) = ( QTRC(k,i,j,I_QI) - QI_t_sat(k,i,j) ) / dt
    enddo
    enddo
    enddo

    call FILE_HISTORY_in( QC_t_sat(:,:,:), 'Pcsat', 'QC production term by satadjust', 'kg/kg/s' )
    call FILE_HISTORY_in( QI_t_sat(:,:,:), 'Pisat', 'QI production term by satadjust', 'kg/kg/s' )

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       EVAPORATE(k,i,j) = max( -QC_t_sat(k,i,j), 0.0_RP ) & ! if negative, condensation
                        * DENS(k,i,j) / (4.0_RP/3.0_RP*PI*DWATR*re_qc**3) ! mass -> number (assuming constant particle radius as re_qc)
    enddo
    enddo
    enddo

    !##### END MP Main #####

    return
  end subroutine ATMOS_PHY_MP_tomita08_adjustment

  !-----------------------------------------------------------------------------
  !> Lin-type cold rain microphysics
  subroutine MP_tomita08( &
         KA, KS, KE, IA, IS, IE, JA, JS, JE, &
         DENS0, PRES0, CCN, &
         dt,                &
         TEMP0, QTRC0,      &
         CPtot0, CVtot0,    &
         RHOE_t             )
    use scale_const, only: &
       UNDEF => CONST_UNDEF, &
       PI    => CONST_PI,    &
       EPS   => CONST_EPS,   &
       Rvap  => CONST_Rvap,  &
       CL    => CONST_CL,    &
       LHV0  => CONST_LHV0,  &
       LHS0  => CONST_LHS0,  &
       LHF0  => CONST_LHF0,  &
       TEM00 => CONST_TEM00, &
       PRE00 => CONST_PRE00, &
       DWATR => CONST_DWATR
    use scale_file_history, only: &
       FILE_HISTORY_query, &
       FILE_HISTORY_put
    use scale_atmos_hydrometeor, only: &
       LHV, &
       LHF, &
       CP_VAPOR, &
       CP_WATER, &
       CP_ICE,   &
       CV_VAPOR, &
       CV_WATER, &
       CV_ICE
    use scale_atmos_saturation, only: &
       SATURATION_dens2qsat_liq => ATMOS_SATURATION_dens2qsat_liq, &
       SATURATION_dens2qsat_ice => ATMOS_SATURATION_dens2qsat_ice
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in) :: DENS0(KA,IA,JA) ! density [km/m3]
    real(RP), intent(in) :: PRES0(KA,IA,JA) ! pressure    [Pa]
    real(RP), intent(in) :: CCN  (KA,IA,JA) ! [1/m3]
    real(DP), intent(in) :: dt

    real(RP), intent(inout) :: TEMP0 (KA,IA,JA) ! temperature [K]
    real(RP), intent(inout) :: QTRC0 (KA,IA,JA,QA_MP)
    real(RP), intent(inout) :: CPtot0(KA,IA,JA)
    real(RP), intent(inout) :: CVtot0(KA,IA,JA)

    real(RP), intent(out) :: RHOE_t(KA,IA,JA)

    real(RP) :: dens(KA)            ! density
    real(RP) :: temp(KA)            ! T [K]
    real(RP) :: cptot, cvtot(KA)
    real(RP) :: qv(KA), qc(KA), qr(KA), qi(KA), qs(KA), qg(KA)
    real(RP) :: qv_t(KA), qc_t(KA), qr_t(KA), qi_t(KA), qs_t(KA), qg_t(KA)
    real(RP) :: e_t, cp_t, cv_t

    real(RP) :: QSATL(KA) ! saturated water vapor for liquid water [kg/kg]
    real(RP) :: QSATI(KA) ! saturated water vapor for ice water    [kg/kg]

    real(RP) :: Sliq(KA)            ! saturated ratio S for liquid water (0-1)
    real(RP) :: Sice(KA)            ! saturated ratio S for ice water    (0-1)

    real(RP) :: Rdens(KA)           ! 1 / density
    real(RP) :: rho_fact(KA)        ! density factor
    real(RP) :: temc(KA)            ! T - T0 [K]

    real(RP) :: N0r(KA), N0s(KA), N0g(KA)

    real(RP) :: RLMDr(KA), RLMDr_2(KA), RLMDr_3(KA)
    real(RP) :: RLMDs, RLMDs_2, RLMDs_3
    real(RP) :: RLMDg(KA), RLMDg_2(KA), RLMDg_3(KA)
    real(RP) :: RLMDr_1br(KA), RLMDr_2br(KA), RLMDr_3br(KA)
    real(RP) :: RLMDs_1bs, RLMDs_2bs, RLMDs_3bs
    real(RP) :: RLMDr_dr(KA), RLMDr_3dr(KA), RLMDr_5dr(KA)
    real(RP) :: RLMDs_ds, RLMDs_3ds, RLMDs_5ds
    real(RP) :: RLMDg_dg(KA), RLMDg_3dg(KA), RLMDg_5dg(KA)
    real(RP) :: RLMDr_7(KA)
    real(RP) :: RLMDr_6dr(KA)

    !---< Roh and Satoh (2014) >---
    real(RP) :: tems, Xs2
    real(RP) :: MOMs_0(KA), MOMs_1(KA), MOMs_2(KA)
    real(RP) :: MOMs_0bs(KA), MOMs_1bs(KA), MOMs_2bs(KA)
    real(RP) :: MOMs_2ds(KA), MOMs_5ds_h(KA), RMOMs_Vt(KA)
    real(RP) :: coef_at(4), coef_bt(4)
    real(RP) :: loga_, b_, nm

    real(RP) :: Vti(KA), Vtr(KA), Vts(KA), Vtg(KA) !< terminal velocity
    real(RP) :: Esi_mod, Egs_mod               !< modified accretion efficiency
    real(RP) :: rhoqc                          !< rho * qc
    real(RP) :: Nc(KA)                         !< Number concentration of cloud water [1/cc]
    real(RP) :: Pracw_orig,  Pracw_kk          !< accretion       term by orig  & k-k scheme
    real(RP) :: Praut_berry, Praut_kk          !< auto-conversion term by berry & k-k scheme
    real(RP) :: Dc                             !< relative variance
    real(RP) :: betai, betas                   !< sticky parameter for auto-conversion
    real(RP) :: Da                             !< thermal diffusion coefficient of air
    real(RP) :: Kd                             !< diffusion coefficient of water vapor in air
    real(RP) :: Nu(KA)                         !< kinematic viscosity of air
    real(RP) :: Glv(KA), Giv(KA), Gil(KA)      !< thermodynamic function
    real(RP) :: ventr, vents, ventg            !< ventilation factor
    real(RP) :: net, fac, fac_sw
    real(RP) :: zerosw, tmp

    !---< Bergeron process >---
    real(RP) :: sw_bergeron(KA)           !< if 0C<T<30C, sw=1
    real(RP) :: a1(KA), a2(KA)
    real(RP) :: ma2(KA) !< 1-a2
    real(RP) :: dt1                   !< time during which the an ice particle of 40um grows to 50um
    real(RP) :: Ni50                  !< number concentration of ice particle of 50um

    !---< Explicit ice generation >---
    real(RP) :: sw, rhoqi, XNi, XMi, Di, Nig, Qig

    logical :: HIST_sw(w_nmax), hist_flag
    real(RP) :: w(KA,w_nmax)

    integer  :: k, i, j, ip
    !---------------------------------------------------------------------------

    call PROF_rapstart('MP_tomita08', 3)

    hist_flag = .false.
    do ip = 1, w_nmax
       call FILE_HISTORY_query( HIST_id(ip), HIST_sw(ip) )
       hist_flag = hist_flag .or. HIST_sw(ip)
    end do

    !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(KA,KS,KE,IS,IE,JS,JE, &
    !$omp        DENS0,TEMP0,PRES0,QTRC0,CCN,CPtot0,CVtot0,dt, &
    !$omp        RHOE_t, &
    !$omp        UNDEF,EPS,PI,PRE00,LHV,LHF,LHF0,CP_VAPOR,CP_WATER,CP_ICE,CV_VAPOR,CV_WATER,CV_ICE,ln10, &
    !$omp        do_couple_aerosol,sw_expice,enable_WDXZ2014,enable_RS2014,enable_KK2000, &
    !$omp        Nc_def,N0r_def,N0s_def,N0g_def, &
    !$omp        Cr,Cs,Cg,Erw,Eri,Eiw,Esw,Esr,Esi,Egw,Egr,Egi,Egs,Ar,As,Ag, &
    !$omp        gamma_sacr,gamma_gacs,gamma_saut,gamma_gaut,beta_saut,beta_gaut,qicrt_saut,qscrt_gaut,mi, &
    !$omp        GAM,GAM_2,GAM_3,GAM_1br,GAM_1bs,GAM_1bsds,GAM_1bg,GAM_1bgdg,GAM_1brdr,GAM_2br,GAM_2bs,GAM_3br,GAM_3bs,GAM_3dr,GAM_3ds,GAM_3dg,GAM_5dr_h,GAM_5ds_h,GAM_5dg_h,GAM_6dr, &
    !$omp        w3d,HIST_sw,hist_flag) &
    !$omp private(dens,temp,cptot,cvtot,qv,qc,qr,qi,qs,qg,qv_t,qc_t,qr_t,qi_t,qs_t,qg_t,e_t,cp_t,cv_t, &
    !$omp         QSATL,QSATI,Sliq,Sice,Rdens,rho_fact,temc,N0r,N0s,N0g, &
    !$omp         RLMDr,RLMDr_2,RLMDr_3,RLMDs,RLMDs_2,RLMDs_3,RLMDg,RLMDg_2,RLMDg_3, &
    !$omp         RLMDr_1br,RLMDr_2br,RLMDr_3br,RLMDs_1bs,RLMDs_2bs,RLMDs_3bs, &
    !$omp         RLMDr_dr,RLMDr_3dr,RLMDr_5dr,RLMDs_ds,RLMDs_3ds,RLMDs_5ds, &
    !$omp         RLMDg_dg,RLMDg_3dg,RLMDg_5dg,RLMDr_7,RLMDr_6dr, &
    !$omp         tems,Xs2,MOMs_0,MOMs_1,MOMs_2,MOMs_0bs,MOMs_1bs,MOMs_2bs,MOMs_2ds,MOMs_5ds_h,RMOMs_Vt, &
    !$omp         coef_at,coef_bt,loga_,b_,nm, &
    !$omp         Vti,Vtr,Vts,Vtg,Esi_mod,Egs_mod,rhoqc,Nc, &
    !$omp         Pracw_orig,Pracw_kk,Praut_berry,Praut_kk,Dc,betai,betas,Da,Kd,Nu, &
    !$omp         Glv,Giv,Gil,ventr,vents,ventg,net,fac,fac_sw,zerosw,tmp, &
    !$omp         sw_bergeron,a1,a2,ma2,dt1,Ni50, &
    !$omp         sw,rhoqi,XNi,XMi,Di,Nig,Qig,w)
!OCL TEMP_PRIVATE(coef_bt,coef_at,w)
    do j = JS, JE
    do i = IS, IE

       if ( do_couple_aerosol ) then
          do k = KS, KE
             Nc(k) = max( CCN(k,i,j)*1.E-6_RP, Nc_def(i,j) ) ! [#/m3]->[#/cc]
          end do
       else
          do k = KS, KE
             Nc(k) = Nc_def(i,j)
          end do
       endif

       ! store to work
       do k = KS, KE
          dens(k) = DENS0(k,i,j)
       end do
       do k = KS, KE
          temp(k) = TEMP0(k,i,j)
       end do

       call SATURATION_dens2qsat_liq( KA, KS, KE, &
                                      temp(:), dens(:), & ! [IN]
                                      QSATL(:)          ) ! [OUT]

       call SATURATION_dens2qsat_ice( KA, KS, KE, &
                                      temp(:), dens(:), & ! [IN]
                                      QSATI(:)          ) ! [OUT]

       ! store to work
       do k = KS, KE
          qv(k) = max( QTRC0(k,i,j,I_QV), 0.0_RP )
       end do
       do k = KS, KE
          qc(k) = max( QTRC0(k,i,j,I_QC), 0.0_RP )
       end do
       do k = KS, KE
          qr(k) = max( QTRC0(k,i,j,I_QR), 0.0_RP )
       end do
       do k = KS, KE
          qi(k) = max( QTRC0(k,i,j,I_QI), 0.0_RP )
       end do
       do k = KS, KE
          qs(k) = max( QTRC0(k,i,j,I_QS), 0.0_RP )
       end do
       do k = KS, KE
          qg(k) = max( QTRC0(k,i,j,I_QG), 0.0_RP )
       end do


       ! saturation ratio S
       do k = KS, KE
          Sliq(k) = qv(k) / max( QSATL(k), EPS )
          Sice(k) = qv(k) / max( QSATI(k), EPS )

          Rdens(k)    = 1.0_RP / dens(k)
          rho_fact(k) = sqrt( dens00 * Rdens(k) )
          temc(k)     = temp(k) - TEM00
       end do

       do k = KS, KE
          w(k,I_delta1) = ( 0.5_RP + sign(0.5_RP, qr(k) - 1.E-4_RP ) )

          w(k,I_delta2) = ( 0.5_RP + sign(0.5_RP, 1.E-4_RP - qr(k) ) ) &
                        * ( 0.5_RP + sign(0.5_RP, 1.E-4_RP - qs(k) ) )

          w(k,I_spsati) = 0.5_RP + sign(0.5_RP, Sice(k) - 1.0_RP )

          w(k,I_iceflg) = 0.5_RP - sign( 0.5_RP, temc(k) ) ! 0: warm, 1: ice
       end do

       do k = KS, KE
          w(k,I_dqv_dt) = qv(k) / dt
          w(k,I_dqc_dt) = qc(k) / dt
          w(k,I_dqr_dt) = qr(k) / dt
          w(k,I_dqi_dt) = qi(k) / dt
          w(k,I_dqs_dt) = qs(k) / dt
          w(k,I_dqg_dt) = qg(k) / dt
       end do

       do k = KS, KE
          sw_bergeron(k) = ( 0.5_RP + sign(0.5_RP, temc(k) + 30.0_RP ) ) &
                         * ( 0.5_RP + sign(0.5_RP, 0.0_RP - temc(k)  ) ) &
                         * ( 1.0_RP - sw_expice                     )
       end do

       ! intercept parameter N0
       if ( enable_WDXZ2014 ) then ! Wainwright et al. (2014)
          do k = KS, KE
             N0r(k) = 1.16E+5_RP * exp( log( max( dens(k)*qr(k)*1000.0_RP, 1.E-2_RP ) )*0.477_RP )
             N0s(k) = 4.58E+9_RP * exp( log( max( dens(k)*qs(k)*1000.0_RP, 1.E-2_RP ) )*0.788_RP )
             N0g(k) = 9.74E+8_RP * exp( log( max( dens(k)*qg(k)*1000.0_RP, 1.E-2_RP ) )*0.816_RP )
          end do
       else
          do k = KS, KE
             N0r(k) = N0r_def ! Marshall and Palmer (1948)
             N0s(k) = N0s_def ! Gunn and Marshall (1958)
             N0g(k) = N0g_def
          end do
       end if

       ! slope parameter lambda (Rain)
       do k = KS, KE
          zerosw = 0.5_RP - sign(0.5_RP, qr(k) - 1.E-12_RP )
          RLMDr    (k) = sqrt(sqrt( dens(k) * qr(k) / ( Ar * N0r(k) * GAM_1br ) + zerosw )) * ( 1.0_RP-zerosw )

          RLMDr_dr (k) = sqrt( RLMDr(k) )       ! **Dr
          RLMDr_2  (k) = RLMDr(k)**2
          RLMDr_3  (k) = RLMDr(k)**3
          RLMDr_7  (k) = RLMDr(k)**7
          RLMDr_1br(k) = RLMDr(k)**4 ! (1+Br)
          RLMDr_2br(k) = RLMDr(k)**5 ! (2+Br)
          RLMDr_3br(k) = RLMDr(k)**6 ! (3+Br)
          RLMDr_3dr(k) = RLMDr(k)**3 * RLMDr_dr(k)
          RLMDr_5dr(k) = RLMDr(k)**5 * RLMDr_dr(k)
          RLMDr_6dr(k) = RLMDr(k)**6 * RLMDr_dr(k)

          w(k,I_RLMDr) = RLMDr(k)
       end do

       ! slope parameter lambda (Snow)
       if ( enable_RS2014 ) then
          !---< modification by Roh and Satoh (2014) >---
          ! bimodal size distribution of snow
          do k = KS, KE
             zerosw = 0.5_RP - sign(0.5_RP, dens(k) * qs(k) - 1.E-12_RP )

             Xs2    = dens(k) * qs(k) / As
             tems       = min( -0.1_RP, temc(k) )
             coef_at(1) = coef_a01 + tems * ( coef_a02 + tems * ( coef_a05 + tems * coef_a09 ) )
             coef_at(2) = coef_a03 + tems * ( coef_a04 + tems *   coef_a07 )
             coef_at(3) = coef_a06 + tems *   coef_a08
             coef_at(4) = coef_a10
             coef_bt(1) = coef_b01 + tems * ( coef_b02 + tems * ( coef_b05 + tems * coef_b09 ) )
             coef_bt(2) = coef_b03 + tems * ( coef_b04 + tems *   coef_b07 )
             coef_bt(3) = coef_b06 + tems *   coef_b08
             coef_bt(4) = coef_b10
             ! 0th moment
             loga_  = coef_at(1)
             b_     = coef_bt(1)
             MOMs_0(k) = exp(ln10*loga_) * exp(log(Xs2+zerosw)*b_) * ( 1.0_RP-zerosw )
             ! 1st moment
             nm = 1.0_RP
             loga_  = coef_at(1) + nm * ( coef_at(2) + nm * ( coef_at(3) + nm * coef_at(4) ) )
             b_  = coef_bt(1) + nm * ( coef_bt(2) + nm * ( coef_bt(3) + nm * coef_bt(4) ) )
             MOMs_1(k) = exp(ln10*loga_) * exp(log(Xs2+zerosw)*b_) * ( 1.0_RP-zerosw )
             ! 2nd moment
             MOMs_2(k) = Xs2
             ! 0 + Bs(=2) moment
             nm = 2.0_RP
             loga_  = coef_at(1) + nm * ( coef_at(2) + nm * ( coef_at(3) + nm * coef_at(4) ) )
             b_  = coef_bt(1) + nm * ( coef_bt(2) + nm * ( coef_bt(3) + nm * coef_bt(4) ) )
             MOMs_0bs(k) = exp(ln10*loga_) * exp(log(Xs2+zerosw)*b_) * ( 1.0_RP-zerosw )
             ! 1 + Bs(=2) moment
             nm = 3.0_RP
             loga_  = coef_at(1) + nm * ( coef_at(2) + nm * ( coef_at(3) + nm * coef_at(4) ) )
             b_  = coef_bt(1) + nm * ( coef_bt(2) + nm * ( coef_bt(3) + nm * coef_bt(4) ) )
             MOMs_1bs(k) = exp(ln10*loga_) * exp(log(Xs2+zerosw)*b_) * ( 1.0_RP-zerosw )
             ! 2 + Bs(=2) moment
             nm = 4.0_RP
             loga_  = coef_at(1) + nm * ( coef_at(2) + nm * ( coef_at(3) + nm * coef_at(4) ) )
             b_  = coef_bt(1) + nm * ( coef_bt(2) + nm * ( coef_bt(3) + nm * coef_bt(4) ) )
             MOMs_2bs(k) = exp(ln10*loga_) * exp(log(Xs2+zerosw)*b_) * ( 1.0_RP-zerosw )
             ! 2 + Ds(=0.25) moment
             nm = 2.25_RP
             loga_  = coef_at(1) + nm * ( coef_at(2) + nm * ( coef_at(3) + nm * coef_at(4) ) )
             b_  = coef_bt(1) + nm * ( coef_bt(2) + nm * ( coef_bt(3) + nm * coef_bt(4) ) )
             MOMs_2ds(k) = exp(ln10*loga_) * exp(log(Xs2+zerosw)*b_) * ( 1.0_RP-zerosw )
             ! ( 3 + Ds(=0.25) ) / 2  moment
             nm = 1.625_RP
             loga_  = coef_at(1) + nm * ( coef_at(2) + nm * ( coef_at(3) + nm * coef_at(4) ) )
             b_  = coef_bt(1) + nm * ( coef_bt(2) + nm * ( coef_bt(3) + nm * coef_bt(4) ) )
             MOMs_5ds_h(k) = exp(ln10*loga_) * exp(log(Xs2+zerosw)*b_) * ( 1.0_RP-zerosw )
             ! Bs(=2) + Ds(=0.25) moment
             nm = 2.25_RP
             loga_  = coef_at(1) + nm * ( coef_at(2) + nm * ( coef_at(3) + nm * coef_at(4) ) )
             b_  = coef_bt(1) + nm * ( coef_bt(2) + nm * ( coef_bt(3) + nm * coef_bt(4) ) )
             RMOMs_Vt(k) = exp(ln10*loga_) * exp(log(Xs2+zerosw)*b_) * ( 1.0_RP-zerosw ) / ( MOMs_0bs(k) + zerosw )

             w(k,I_RLMDs) = UNDEF
          end do
       else
          do k = KS, KE
             zerosw = 0.5_RP - sign(0.5_RP, dens(k) * qs(k) - 1.E-12_RP )

             RLMDs  = sqrt(sqrt( dens(k) * qs(k) / ( As * N0s(k) * GAM_1bs ) + zerosw )) * ( 1.0_RP-zerosw )
             RLMDs_ds  = sqrt( sqrt(RLMDs) ) ! **Ds
             RLMDs_2   = RLMDs**2
             RLMDs_3   = RLMDs_2 * RLMDs
             RLMDs_1bs = RLMDs_2 * RLMDs_2 ! (1+Bs)
             RLMDs_2bs = RLMDs_3 * RLMDs_2 ! (2+Bs)
             RLMDs_3bs = RLMDs_3 * RLMDs_3 ! (3+Bs)
             RLMDs_3ds = RLMDs_3 * RLMDs_ds
             RLMDs_5ds = RLMDs_2 * RLMDs_3ds

             MOMs_0    (k) = N0s(k) * GAM       * RLMDs           ! Ns * 0th moment
             MOMs_1    (k) = N0s(k) * GAM_2     * RLMDs_2         ! Ns * 1st moment
             MOMs_2    (k) = N0s(k) * GAM_3     * RLMDs_3         ! Ns * 2nd moment
             MOMs_0bs  (k) = N0s(k) * GAM_1bs   * RLMDs_1bs       ! Ns * 0+bs
             MOMs_1bs  (k) = N0s(k) * GAM_2bs   * RLMDs_2bs       ! Ns * 1+bs
             MOMs_2bs  (k) = N0s(k) * GAM_3bs   * RLMDs_3bs       ! Ns * 2+bs
             MOMs_2ds  (k) = N0s(k) * GAM_3ds   * RLMDs_3ds       ! Ns * 2+ds
             MOMs_5ds_h(k) = N0s(k) * GAM_5ds_h * sqrt(RLMDs_5ds) ! Ns * (5+ds)/2
             RMOMs_Vt  (k) = GAM_1bsds / GAM_1bs * RLMDs_ds

             w(k,I_RLMDs) = RLMDs
          end do
       end if

       ! slope parameter lambda (Graupel)
       do k = KS, KE
          zerosw = 0.5_RP - sign(0.5_RP, qg(k) - 1.E-12_RP )
          RLMDg(k)  = sqrt(sqrt( dens(k) * qg(k) / ( Ag * N0g(k) * GAM_1bg ) + zerosw )) * ( 1.0_RP-zerosw )

          RLMDg_dg (k) = sqrt( RLMDg(k) )       ! **Dg
          RLMDg_2  (k) = RLMDg(k)**2
          RLMDg_3  (k) = RLMDg(k) * RLMDg_2(k)
          RLMDg_3dg(k) = RLMDg_3(k) * RLMDg_dg(k)
          RLMDg_5dg(k) = RLMDg_2(k) * RLMDg_3dg(k)

          w(k,I_RLMDg) = RLMDg(k)
       end do

       !---< terminal velocity >
       do k = KS, KE
          zerosw = 0.5_RP - sign(0.5_RP, qi(k) - 1.E-8_RP )
          Vti(k) = -3.29_RP * exp( log( dens(k)*qi(k)+zerosw )*0.16_RP ) * ( 1.0_RP-zerosw )
          Vtr(k) = -Cr * rho_fact(k) * GAM_1brdr / GAM_1br * RLMDr_dr(k)
          Vts(k) = -Cs * rho_fact(k) * RMOMs_Vt(k)
          Vtg(k) = -Cg * rho_fact(k) * GAM_1bgdg / GAM_1bg * RLMDg_dg(k)
       end do

       !---< Nucleation >---
       ! [Pigen] ice nucleation
       do k = KS, KE
          Nig = max( exp(-0.1_RP*temc(k)), 1.0_RP ) * 1000.0_RP
          Qig = 4.92E-11_RP * exp(log(Nig)*1.33_RP) * Rdens(k)

          w(k,I_Pigen) = max( min( Qig-qi(k), qv(k)-QSATI(k) ), 0.0_RP ) / dt
       end do

       !---< Accretion >---

       ! [Pracw] accretion rate of cloud water by rain
       if ( enable_KK2000 ) then
          do k = KS, KE
             zerosw     = 0.5_RP - sign(0.5_RP, qc(k)*qr(k) - 1.E-12_RP )
             Pracw_kk   = 67.0_RP * exp( log( qc(k)*qr(k)+zerosw )*1.15_RP ) * ( 1.0_RP-zerosw ) ! eq.(33) in KK(2000)
             w(k,I_Pracw) = Pracw_kk
          end do
       else
          do k = KS, KE
             Pracw_orig = qc(k) * 0.25_RP * PI * Erw * N0r(k) * Cr * GAM_3dr * RLMDr_3dr(k) * rho_fact(k)
             w(k,I_Pracw) = Pracw_orig
          end do
       end if


       do k = KS, KE
          ! [Psacw] accretion rate of cloud water by snow
          w(k,I_Psacw) = qc(k) * 0.25_RP * PI * Esw          * Cs * MOMs_2ds(k)            * rho_fact(k)
          ! [Pgacw] accretion rate of cloud water by graupel
          w(k,I_Pgacw) = qc(k) * 0.25_RP * PI * Egw * N0g(k) * Cg * GAM_3dg * RLMDg_3dg(k) * rho_fact(k)
       end do

       do k = KS, KE
          Esi_mod = min( Esi, Esi * exp( gamma_sacr * temc(k) ) )
          ! [Praci] accretion rate of cloud ice by rain
          w(k,I_Praci) = qi(k) * 0.25_RP * PI * Eri * N0r(k) * Cr * GAM_3dr * RLMDr_3dr(k) * rho_fact(k)
          ! [Psaci] accretion rate of cloud ice by snow
          w(k,I_Psaci) = qi(k) * 0.25_RP * PI * Esi_mod      * Cs * MOMs_2ds(k)            * rho_fact(k)
          ! [Pgaci] accretion rate of cloud ice by grupel
          w(k,I_Pgaci) = qi(k) * 0.25_RP * PI * Egi * N0g(k) * Cg * GAM_3dg * RLMDg_3dg(k) * rho_fact(k)
          ! [Piacr] accretion rate of rain by cloud ice
          w(k,I_Piacr) = qi(k) * Ar / mi * 0.25_RP * PI * Eri * N0r(k) * Cr * GAM_6dr * RLMDr_6dr(k) * rho_fact(k)
       end do

       do k = KS, KE
          ! [Psacr] accretion rate of rain by snow
          w(k,I_Psacr) = Ar * 0.25_RP * PI * Rdens(k) * Esr * N0r(k)          * abs(Vtr(k)-Vts(k)) &
                       * (          GAM_1br * RLMDr_1br(k) * MOMs_2(k)          &
                         + 2.0_RP * GAM_2br * RLMDr_2br(k) * MOMs_1(k)          &
                         +          GAM_3br * RLMDr_3br(k) * MOMs_0(k)          )

          ! [Pgacr] accretion rate of rain by graupel
          w(k,I_Pgacr) = Ar * 0.25_RP * PI * Rdens(k) * Egr * N0g(k) * N0r(k) * abs(Vtg(k)-Vtr(k)) &
                       * (          GAM_1br * RLMDr_1br(k) * GAM_3 * RLMDg_3(k) &
                         + 2.0_RP * GAM_2br * RLMDr_2br(k) * GAM_2 * RLMDg_2(k) &
                         +          GAM_3br * RLMDr_3br(k) * GAM   * RLMDg  (k) )

          ! [Pracs] accretion rate of snow by rain
          w(k,I_Pracs) = As * 0.25_RP * PI * Rdens(k) * Esr       *  N0r(k)   * abs(Vtr(k)-Vts(k)) &
                       * (          MOMs_0bs(k)            * GAM_3 * RLMDr_3(k) &
                         + 2.0_RP * MOMs_1bs(k)            * GAM_2 * RLMDr_2(k) &
                         +          MOMs_2bs(k)            * GAM   * RLMDr  (k) )

          ! [Pgacs] accretion rate of snow by graupel
          Egs_mod = min( Egs, Egs * exp( gamma_gacs * temc(k) ) )
          w(k,I_Pgacs) = As * 0.25_RP * PI * Rdens(k) * Egs_mod   * N0g(k)   * abs(Vtg(k)-Vts(k)) &
                       * (          MOMs_0bs(k)            * GAM_3 * RLMDg_3(k) &
                         + 2.0_RP * MOMs_1bs(k)            * GAM_2 * RLMDg_2(k) &
                         +          MOMs_2bs(k)            * GAM   * RLMDg  (k) )
       end do

       !---< Auto-conversion >---
       ! [Praut] auto-conversion rate from cloud water to rain
       if ( enable_KK2000 ) then
          do k = KS, KE
             zerosw      = 0.5_RP - sign(0.5_RP, qc(k) - 1.E-12_RP )
             Praut_kk    = 1350.0_RP                                              &
                         * exp( log( qc(k)+zerosw )*2.47_RP ) * ( 1.0_RP-zerosw ) &
                         * exp( log( Nc(k) )*(-1.79_RP) )                    ! eq.(29) in KK(2000)
             w(k,I_Praut) = Praut_kk
          end do
       else
          do k = KS, KE
             rhoqc = dens(k) * qc(k) * 1000.0_RP ! [g/m3]
             Dc    = 0.146_RP - 5.964E-2_RP * log( Nc(k) / 2000.0_RP )
             Praut_berry = Rdens(k) * 1.67E-5_RP * rhoqc * rhoqc / ( 5.0_RP + 3.66E-2_RP * Nc(k) / ( Dc * rhoqc + EPS ) )
             w(k,I_Praut) = Praut_berry
          end do
       end if

       do k = KS, KE
          ! [Psaut] auto-conversion rate from cloud ice to snow
          betai = min( beta_saut, beta_saut * exp( gamma_saut * temc(k) ) )
          w(k,I_Psaut) = max( betai*(qi(k)-qicrt_saut), 0.0_RP )
          ! [Pgaut] auto-conversion rate from snow to graupel
          betas = min( beta_gaut, beta_gaut * exp( gamma_gaut * temc(k) ) )
          w(k,I_Pgaut) = max( betas*(qs(k)-qscrt_gaut), 0.0_RP )
       end do

       !---< Evaporation, Sublimation, Melting, and Freezing >---
       do k = KS, KE
          Da    = ( Da0 + dDa_dT * temc(k) )
          Kd    = ( Dw0 + dDw_dT * temc(k) ) * PRE00 / PRES0(k,i,j)
          NU(k) = ( mu0 + dmu_dT * temc(k) ) * Rdens(k)

          Glv(k) = 1.0_RP / ( LHV0/(Da*temp(k)) * ( LHV0/(Rvap*temp(k)) - 1.0_RP ) + 1.0_RP/(Kd*dens(k)*QSATL(k)) )
          Giv(k) = 1.0_RP / ( LHS0/(Da*temp(k)) * ( LHS0/(Rvap*temp(k)) - 1.0_RP ) + 1.0_RP/(Kd*dens(k)*QSATI(k)) )
          Gil(k) = ( Da * temc(k) ) / LHF0
       end do

       ! [Prevp] evaporation rate of rain
       do k = KS, KE
          ventr = f1r * GAM_2 * RLMDr_2(k) + f2r * sqrt( Cr * rho_fact(k) / NU(k) * RLMDr_5dr(k) ) * GAM_5dr_h
          w(k,I_Prevp) = 2.0_RP * PI * Rdens(k) * N0r(k) * ( 1.0_RP-min(Sliq(k),1.0_RP) ) * Glv(k) * ventr
       end do

       ! [Pidep,Pisub] deposition/sublimation rate for ice
       do k = KS, KE
          rhoqi = max(dens(k)*qi(k), EPS)
          XNi   = min( max( 5.38E+7_RP * exp( log(rhoqi)*0.75_RP ), 1.E+3_RP ), 1.E+6_RP )
          XMi   = rhoqi / XNi
          Di    = min( Di_a * sqrt(XMi), Di_max )
          tmp = 4.0_RP * Di * XNi * Rdens(k) * ( Sice(k)-1.0_RP ) * Giv(k)
          w(k,I_Pidep) = (        w(k,I_spsati) ) * ( tmp) ! Sice > 1
          w(k,I_Pisub) = ( 1.0_RP-w(k,I_spsati) ) * (-tmp) ! Sice < 1
       end do

       ! [Pihom] homogenious freezing at T < -40C
       do k = KS, KE
          sw = ( 0.5_RP - sign(0.5_RP, temc(k) + 40.0_RP ) ) ! if T < -40C, sw=1
          w(k,I_Pihom) = sw * qc(k) / dt
       end do

       ! [Pihtr] heteroginous freezing at -40C < T < 0C
       do k = KS, KE
          sw = ( 0.5_RP + sign(0.5_RP, temc(k) + 40.0_RP ) ) &
             * ( 0.5_RP - sign(0.5_RP, temc(k)           ) ) ! if -40C < T < 0C, sw=1
          w(k,I_Pihtr) = sw * ( dens(k) / DWATR * qc(k)**2 / ( Nc_ihtr * 1.E+6_RP ) ) &
                       * B_frz * ( exp(-A_frz*temc(k)) - 1.0_RP )
       end do

       ! [Pimlt] ice melting at T > 0C
       do k = KS, KE
          sw = ( 0.5_RP + sign(0.5_RP, temc(k)           ) ) ! if T > 0C, sw=1
          w(k,I_Pimlt) = sw * qi(k) / dt
       end do

       do k = KS, KE
          ! [Psdep,Pssub] deposition/sublimation rate for snow
          vents = f1s * MOMs_1(k)          + f2s * sqrt( Cs * rho_fact(k) / NU(k)             ) * MOMs_5ds_h(k)
          tmp = 2.0_RP * PI * Rdens(k) *       ( Sice(k)-1.0_RP ) * Giv(k) * vents
          w(k,I_Psdep) = (        w(k,I_spsati) ) * ( tmp) ! Sice > 1
          w(k,I_Pssub) = ( 1.0_RP-w(k,I_spsati) ) * (-tmp) ! Sice < 1
          ! [Psmlt] melting rate of snow
          w(k,I_Psmlt) = 2.0_RP * PI * Rdens(k) *       Gil(k) * vents &
                       + CL * temc(k) / LHF0 * ( w(k,I_Psacw) + w(k,I_Psacr) )
          w(k,I_Psmlt) = max( w(k,I_Psmlt), 0.0_RP )
       end do

       do k = KS, KE
          ! [Pgdep/pgsub] deposition/sublimation rate for graupel
          ventg = f1g * GAM_2 * RLMDg_2(k) + f2g * sqrt( Cg * rho_fact(k) / NU(k) * RLMDg_5dg(k) ) * GAM_5dg_h
          tmp = 2.0_RP * PI * Rdens(k) * N0g(k) * ( Sice(k)-1.0_RP ) * Giv(k) * ventg
          w(k,I_Pgdep) = (        w(k,I_spsati) ) * ( tmp) ! Sice > 1
          w(k,I_Pgsub) = ( 1.0_RP-w(k,I_spsati) ) * (-tmp) ! Sice < 1
          ! [Pgmlt] melting rate of graupel
          w(k,I_Pgmlt) = 2.0_RP * PI * Rdens(k) * N0g(k) * Gil(k) * ventg &
                       + CL * temc(k) / LHF0 * ( w(k,I_Pgacw) + w(k,I_Pgacr) )
          w(k,I_Pgmlt) = max( w(k,I_Pgmlt), 0.0_RP )
       end do

       ! [Pgfrz] freezing rate of graupel
       do k = KS, KE
          w(k,I_Pgfrz) = 2.0_RP * PI * Rdens(k) * N0r(k) * 60.0_RP * B_frz * Ar * ( exp(-A_frz*temc(k)) - 1.0_RP ) * RLMDr_7(k)
       end do

       ! [Psfw,Psfi] ( Bergeron process ) growth rate of snow by Bergeron process from cloud water/ice
       call MP_tomita08_BergeronParam( KA, KS, KE, &
                                       temp(:),             & ! [IN]
                                       a1(:), a2(:), ma2(:) ) ! [OUT]
       do k = KS, KE
          dt1  = ( exp( log(mi50)*ma2(k) ) &
                 - exp( log(mi40)*ma2(k) ) ) / ( a1(k) * ma2(k) )
          Ni50 = qi(k) * dt / ( mi50 * dt1 )
          w(k,I_Psfw ) = Ni50 * ( a1(k) * exp( log(mi50)*a2(k) )                 &
                                + PI * Eiw * dens(k) * qc(k) * Ri50*Ri50 * vti50 )
          w(k,I_Psfi ) = qi(k) / dt1
       end do

       !---< limiter >---
       do k = KS, KE
          w(k,I_Pigen) = min( w(k,I_Pigen), w(k,I_dqv_dt) ) * (        w(k,I_iceflg) ) * sw_expice
          w(k,I_Pidep) = min( w(k,I_Pidep), w(k,I_dqv_dt) ) * (        w(k,I_iceflg) ) * sw_expice
          w(k,I_Psdep) = min( w(k,I_Psdep), w(k,I_dqv_dt) ) * (        w(k,I_iceflg) )
          w(k,I_Pgdep) = min( w(k,I_Pgdep), w(k,I_dqv_dt) ) * (        w(k,I_iceflg) )
       end do

       do k = KS, KE
          w(k,I_Pracw) = w(k,I_Pracw)                          &
                       + w(k,I_Psacw) * ( 1.0_RP-w(k,I_iceflg) ) & ! c->r by s
                       + w(k,I_Pgacw) * ( 1.0_RP-w(k,I_iceflg) )   ! c->r by g
       end do

       do k = KS, KE
          w(k,I_Praut) = min( w(k,I_Praut), w(k,I_dqc_dt) )
          w(k,I_Pracw) = min( w(k,I_Pracw), w(k,I_dqc_dt) )
          w(k,I_Pihom) = min( w(k,I_Pihom), w(k,I_dqc_dt) ) * (        w(k,I_iceflg) ) * sw_expice
          w(k,I_Pihtr) = min( w(k,I_Pihtr), w(k,I_dqc_dt) ) * (        w(k,I_iceflg) ) * sw_expice
          w(k,I_Psacw) = min( w(k,I_Psacw), w(k,I_dqc_dt) ) * (        w(k,I_iceflg) )
          w(k,I_Psfw ) = min( w(k,I_Psfw ), w(k,I_dqc_dt) ) * (        w(k,I_iceflg) ) * sw_bergeron(k)
          w(k,I_Pgacw) = min( w(k,I_Pgacw), w(k,I_dqc_dt) ) * (        w(k,I_iceflg) )
       end do

       do k = KS, KE
          w(k,I_Prevp) = min( w(k,I_Prevp), w(k,I_dqr_dt) )
          w(k,I_Piacr) = min( w(k,I_Piacr), w(k,I_dqr_dt) ) * (        w(k,I_iceflg) )
          w(k,I_Psacr) = min( w(k,I_Psacr), w(k,I_dqr_dt) ) * (        w(k,I_iceflg) )
          w(k,I_Pgacr) = min( w(k,I_Pgacr), w(k,I_dqr_dt) ) * (        w(k,I_iceflg) )
          w(k,I_Pgfrz) = min( w(k,I_Pgfrz), w(k,I_dqr_dt) ) * (        w(k,I_iceflg) )
       end do

       do k = KS, KE
          w(k,I_Pisub) = min( w(k,I_Pisub), w(k,I_dqi_dt) ) * (        w(k,I_iceflg) ) * sw_expice
          w(k,I_Pimlt) = min( w(k,I_Pimlt), w(k,I_dqi_dt) ) * ( 1.0_RP-w(k,I_iceflg) ) * sw_expice
          w(k,I_Psaut) = min( w(k,I_Psaut), w(k,I_dqi_dt) ) * (        w(k,I_iceflg) )
          w(k,I_Praci) = min( w(k,I_Praci), w(k,I_dqi_dt) ) * (        w(k,I_iceflg) )
          w(k,I_Psaci) = min( w(k,I_Psaci), w(k,I_dqi_dt) ) * (        w(k,I_iceflg) )
          w(k,I_Psfi ) = min( w(k,I_Psfi ), w(k,I_dqi_dt) ) * (        w(k,I_iceflg) ) * sw_bergeron(k)
          w(k,I_Pgaci) = min( w(k,I_Pgaci), w(k,I_dqi_dt) ) * (        w(k,I_iceflg) )
       end do

       do k = KS, KE
          w(k,I_Pssub) = min( w(k,I_Pssub), w(k,I_dqs_dt) ) * (        w(k,I_iceflg) )
          w(k,I_Psmlt) = min( w(k,I_Psmlt), w(k,I_dqs_dt) ) * ( 1.0_RP-w(k,I_iceflg) )
          w(k,I_Pgaut) = min( w(k,I_Pgaut), w(k,I_dqs_dt) ) * (        w(k,I_iceflg) )
          w(k,I_Pracs) = min( w(k,I_Pracs), w(k,I_dqs_dt) ) * (        w(k,I_iceflg) )
          w(k,I_Pgacs) = min( w(k,I_Pgacs), w(k,I_dqs_dt) )
       end do

       do k = KS, KE
          w(k,I_Pgsub) = min( w(k,I_Pgsub), w(k,I_dqg_dt) ) * (        w(k,I_iceflg) )
          w(k,I_Pgmlt) = min( w(k,I_Pgmlt), w(k,I_dqg_dt) ) * ( 1.0_RP-w(k,I_iceflg) )
       end do

       do k = KS, KE
          w(k,I_Piacr_s) = ( 1.0_RP - w(k,I_delta1) ) * w(k,I_Piacr)
          w(k,I_Piacr_g) = (          w(k,I_delta1) ) * w(k,I_Piacr)
          w(k,I_Praci_s) = ( 1.0_RP - w(k,I_delta1) ) * w(k,I_Praci)
          w(k,I_Praci_g) = (          w(k,I_delta1) ) * w(k,I_Praci)
          w(k,I_Psacr_s) = (          w(k,I_delta2) ) * w(k,I_Psacr)
          w(k,I_Psacr_g) = ( 1.0_RP - w(k,I_delta2) ) * w(k,I_Psacr)
          w(k,I_Pracs  ) = ( 1.0_RP - w(k,I_delta2) ) * w(k,I_Pracs)
       end do

       ! [QC]
       do k = KS, KE
          net = &
              + w(k,I_Pimlt  ) & ! [prod] i->c
              - w(k,I_Praut  ) & ! [loss] c->r
              - w(k,I_Pracw  ) & ! [loss] c->r
              - w(k,I_Pihom  ) & ! [loss] c->i
              - w(k,I_Pihtr  ) & ! [loss] c->i
              - w(k,I_Psacw  ) & ! [loss] c->s
              - w(k,I_Psfw   ) & ! [loss] c->s
              - w(k,I_Pgacw  )   ! [loss] c->g

          fac_sw = 0.5_RP + sign( 0.5_RP, net+EPS ) ! if production > loss , fac_sw=1
          fac    = (          fac_sw ) &
                 + ( 1.0_RP - fac_sw ) * min( -w(k,I_dqc_dt)/(net-fac_sw), 1.0_RP ) ! loss limiter

          w(k,I_Pimlt  ) = w(k,I_Pimlt  ) * fac
          w(k,I_Praut  ) = w(k,I_Praut  ) * fac
          w(k,I_Pracw  ) = w(k,I_Pracw  ) * fac
          w(k,I_Pihom  ) = w(k,I_Pihom  ) * fac
          w(k,I_Pihtr  ) = w(k,I_Pihtr  ) * fac
          w(k,I_Psacw  ) = w(k,I_Psacw  ) * fac
          w(k,I_Psfw   ) = w(k,I_Psfw   ) * fac
          w(k,I_Pgacw  ) = w(k,I_Pgacw  ) * fac
       end do

       ! [QI]
       do k = KS, KE
          net = &
              + w(k,I_Pigen  ) & ! [prod] v->i
              + w(k,I_Pidep  ) & ! [prod] v->i
              + w(k,I_Pihom  ) & ! [prod] c->i
              + w(k,I_Pihtr  ) & ! [prod] c->i
              - w(k,I_Pisub  ) & ! [loss] i->v
              - w(k,I_Pimlt  ) & ! [loss] i->c
              - w(k,I_Psaut  ) & ! [loss] i->s
              - w(k,I_Praci_s) & ! [loss] i->s
              - w(k,I_Psaci  ) & ! [loss] i->s
              - w(k,I_Psfi   ) & ! [loss] i->s
              - w(k,I_Praci_g) & ! [loss] i->g
              - w(k,I_Pgaci  )   ! [loss] i->g

          fac_sw = 0.5_RP + sign( 0.5_RP, net+EPS ) ! if production > loss , fac_sw=1
          fac    = (          fac_sw ) &
                 + ( 1.0_RP - fac_sw ) * min( -w(k,I_dqi_dt)/(net-fac_sw), 1.0_RP ) ! loss limiter

          w(k,I_Pigen  ) = w(k,I_Pigen  ) * fac
          w(k,I_Pidep  ) = w(k,I_Pidep  ) * fac
          w(k,I_Pihom  ) = w(k,I_Pihom  ) * fac
          w(k,I_Pihtr  ) = w(k,I_Pihtr  ) * fac
          w(k,I_Pisub  ) = w(k,I_Pisub  ) * fac
          w(k,I_Pimlt  ) = w(k,I_Pimlt  ) * fac
          w(k,I_Psaut  ) = w(k,I_Psaut  ) * fac
          w(k,I_Praci_s) = w(k,I_Praci_s) * fac
          w(k,I_Psaci  ) = w(k,I_Psaci  ) * fac
          w(k,I_Psfi   ) = w(k,I_Psfi   ) * fac
          w(k,I_Praci_g) = w(k,I_Praci_g) * fac
          w(k,I_Pgaci  ) = w(k,I_Pgaci  ) * fac
       end do

       ! [QR]
       do k = KS, KE
          net = &
              + w(k,I_Praut  ) & ! [prod] c->r
              + w(k,I_Pracw  ) & ! [prod] c->r
              + w(k,I_Psmlt  ) & ! [prod] s->r
              + w(k,I_Pgmlt  ) & ! [prod] g->r
              - w(k,I_Prevp  ) & ! [loss] r->v
              - w(k,I_Piacr_s) & ! [loss] r->s
              - w(k,I_Psacr_s) & ! [loss] r->s
              - w(k,I_Piacr_g) & ! [loss] r->g
              - w(k,I_Psacr_g) & ! [loss] r->g
              - w(k,I_Pgacr  ) & ! [loss] r->g
              - w(k,I_Pgfrz  )   ! [loss] r->g

          fac_sw = 0.5_RP + sign( 0.5_RP, net+EPS ) ! if production > loss , fac_sw=1
          fac    = (          fac_sw ) &
                 + ( 1.0_RP - fac_sw ) * min( -w(k,I_dqr_dt)/(net-fac_sw), 1.0_RP ) ! loss limiter

          w(k,I_Praut  ) = w(k,I_Praut  ) * fac
          w(k,I_Pracw  ) = w(k,I_Pracw  ) * fac
          w(k,I_Psmlt  ) = w(k,I_Psmlt  ) * fac
          w(k,I_Pgmlt  ) = w(k,I_Pgmlt  ) * fac
          w(k,I_Prevp  ) = w(k,I_Prevp  ) * fac
          w(k,I_Piacr_s) = w(k,I_Piacr_s) * fac
          w(k,I_Psacr_s) = w(k,I_Psacr_s) * fac
          w(k,I_Piacr_g) = w(k,I_Piacr_g) * fac
          w(k,I_Psacr_g) = w(k,I_Psacr_g) * fac
          w(k,I_Pgacr  ) = w(k,I_Pgacr  ) * fac
          w(k,I_Pgfrz  ) = w(k,I_Pgfrz  ) * fac
       end do

       ! [QV]
       do k = KS, KE
          net = &
              + w(k,I_Prevp  ) & ! [prod] r->v
              + w(k,I_Pisub  ) & ! [prod] i->v
              + w(k,I_Pssub  ) & ! [prod] s->v
              + w(k,I_Pgsub  ) & ! [prod] g->v
              - w(k,I_Pigen  ) & ! [loss] v->i
              - w(k,I_Pidep  ) & ! [loss] v->i
              - w(k,I_Psdep  ) & ! [loss] v->s
              - w(k,I_Pgdep  )   ! [loss] v->g

          fac_sw = 0.5_RP + sign( 0.5_RP, net+EPS ) ! if production > loss , fac_sw=1
          fac    = (          fac_sw ) &
                 + ( 1.0_RP - fac_sw ) * min( -w(k,I_dqv_dt)/(net-fac_sw), 1.0_RP ) ! loss limiter

          w(k,I_Prevp  ) = w(k,I_Prevp  ) * fac
          w(k,I_Pisub  ) = w(k,I_Pisub  ) * fac
          w(k,I_Pssub  ) = w(k,I_Pssub  ) * fac
          w(k,I_Pgsub  ) = w(k,I_Pgsub  ) * fac
          w(k,I_Pigen  ) = w(k,I_Pigen  ) * fac
          w(k,I_Pidep  ) = w(k,I_Pidep  ) * fac
          w(k,I_Psdep  ) = w(k,I_Psdep  ) * fac
          w(k,I_Pgdep  ) = w(k,I_Pgdep  ) * fac
       end do

       ! [QS]
       do k = KS, KE
          net = &
              + w(k,I_Psdep  ) & ! [prod] v->s
              + w(k,I_Psacw  ) & ! [prod] c->s
              + w(k,I_Psfw   ) & ! [prod] c->s
              + w(k,I_Piacr_s) & ! [prod] r->s
              + w(k,I_Psacr_s) & ! [prod] r->s
              + w(k,I_Psaut  ) & ! [prod] i->s
              + w(k,I_Praci_s) & ! [prod] i->s
              + w(k,I_Psaci  ) & ! [prod] i->s
              + w(k,I_Psfi   ) & ! [prod] i->s
              - w(k,I_Pssub  ) & ! [loss] s->v
              - w(k,I_Psmlt  ) & ! [loss] s->r
              - w(k,I_Pgaut  ) & ! [loss] s->g
              - w(k,I_Pracs  ) & ! [loss] s->g
              - w(k,I_Pgacs  )   ! [loss] s->g

          fac_sw = 0.5_RP + sign( 0.5_RP, net+EPS ) ! if production > loss , fac_sw=1
          fac    = (          fac_sw ) &
                 + ( 1.0_RP - fac_sw ) * min( -w(k,I_dqs_dt)/(net-fac_sw), 1.0_RP ) ! loss limiter

          w(k,I_Psdep  ) = w(k,I_Psdep  ) * fac
          w(k,I_Psacw  ) = w(k,I_Psacw  ) * fac
          w(k,I_Psfw   ) = w(k,I_Psfw   ) * fac
          w(k,I_Piacr_s) = w(k,I_Piacr_s) * fac
          w(k,I_Psacr_s) = w(k,I_Psacr_s) * fac
          w(k,I_Psaut  ) = w(k,I_Psaut  ) * fac
          w(k,I_Praci_s) = w(k,I_Praci_s) * fac
          w(k,I_Psaci  ) = w(k,I_Psaci  ) * fac
          w(k,I_Psfi   ) = w(k,I_Psfi   ) * fac
          w(k,I_Pssub  ) = w(k,I_Pssub  ) * fac
          w(k,I_Psmlt  ) = w(k,I_Psmlt  ) * fac
          w(k,I_Pgaut  ) = w(k,I_Pgaut  ) * fac
          w(k,I_Pracs  ) = w(k,I_Pracs  ) * fac
          w(k,I_Pgacs  ) = w(k,I_Pgacs  ) * fac
       end do

       ! [QG]
       do k = KS, KE
          net = &
              + w(k,I_Pgdep  ) & ! [prod] v->g
              + w(k,I_Pgacw  ) & ! [prod] c->g
              + w(k,I_Piacr_g) & ! [prod] r->g
              + w(k,I_Psacr_g) & ! [prod] r->g
              + w(k,I_Pgacr  ) & ! [prod] r->g
              + w(k,I_Pgfrz  ) & ! [prod] r->g
              + w(k,I_Praci_g) & ! [prod] i->g
              + w(k,I_Pgaci  ) & ! [prod] i->g
              + w(k,I_Pgaut  ) & ! [prod] s->g
              + w(k,I_Pracs  ) & ! [prod] s->g
              + w(k,I_Pgacs  ) & ! [prod] s->g
              - w(k,I_Pgsub  ) & ! [loss] g->v
              - w(k,I_Pgmlt  )   ! [loss] g->r

          fac_sw = 0.5_RP + sign( 0.5_RP, net+EPS ) ! if production > loss , fac_sw=1
          fac    = (          fac_sw ) &
                 + ( 1.0_RP - fac_sw ) * min( -w(k,I_dqg_dt)/(net-fac_sw), 1.0_RP ) ! loss limiter

          w(k,I_Pgdep  ) = w(k,I_Pgdep  ) * fac
          w(k,I_Pgacw  ) = w(k,I_Pgacw  ) * fac
          w(k,I_Piacr_g) = w(k,I_Piacr_g) * fac
          w(k,I_Psacr_g) = w(k,I_Psacr_g) * fac
          w(k,I_Pgacr  ) = w(k,I_Pgacr  ) * fac
          w(k,I_Pgfrz  ) = w(k,I_Pgfrz  ) * fac
          w(k,I_Praci_g) = w(k,I_Praci_g) * fac
          w(k,I_Pgaci  ) = w(k,I_Pgaci  ) * fac
          w(k,I_Pgaut  ) = w(k,I_Pgaut  ) * fac
          w(k,I_Pracs  ) = w(k,I_Pracs  ) * fac
          w(k,I_Pgacs  ) = w(k,I_Pgacs  ) * fac
          w(k,I_Pgsub  ) = w(k,I_Pgsub  ) * fac
          w(k,I_Pgmlt  ) = w(k,I_Pgmlt  ) * fac
       end do


       do k = KS, KE
          qc_t(k) = + w(k,I_Pimlt  ) & ! [prod] i->c
                   - w(k,I_Praut  ) & ! [loss] c->r
                   - w(k,I_Pracw  ) & ! [loss] c->r
                   - w(k,I_Pihom  ) & ! [loss] c->i
                   - w(k,I_Pihtr  ) & ! [loss] c->i
                   - w(k,I_Psacw  ) & ! [loss] c->s
                   - w(k,I_Psfw   ) & ! [loss] c->s
                   - w(k,I_Pgacw  )   ! [loss] c->g
          qc_t(k) = max( qc_t(k), -w(k,I_dqc_dt) )
       end do

       do k = KS, KE
          qr_t(k) = + w(k,I_Praut  ) & ! [prod] c->r
                    + w(k,I_Pracw  ) & ! [prod] c->r
                    + w(k,I_Psmlt  ) & ! [prod] s->r
                    + w(k,I_Pgmlt  ) & ! [prod] g->r
                    - w(k,I_Prevp  ) & ! [loss] r->v
                    - w(k,I_Piacr_s) & ! [loss] r->s
                    - w(k,I_Psacr_s) & ! [loss] r->s
                    - w(k,I_Piacr_g) & ! [loss] r->g
                    - w(k,I_Psacr_g) & ! [loss] r->g
                    - w(k,I_Pgacr  ) & ! [loss] r->g
                    - w(k,I_Pgfrz  )   ! [loss] r->g

          qr_t(k) = max( qr_t(k), -w(k,I_dqr_dt) )
       end do

       do k = KS, KE
          qi_t(k) = + w(k,I_Pigen  ) & ! [prod] v->i
                    + w(k,I_Pidep  ) & ! [prod] v->i
                    + w(k,I_Pihom  ) & ! [prod] c->i
                    + w(k,I_Pihtr  ) & ! [prod] c->i
                    - w(k,I_Pisub  ) & ! [loss] i->v
                    - w(k,I_Pimlt  ) & ! [loss] i->c
                    - w(k,I_Psaut  ) & ! [loss] i->s
                    - w(k,I_Praci_s) & ! [loss] i->s
                    - w(k,I_Psaci  ) & ! [loss] i->s
                    - w(k,I_Psfi   ) & ! [loss] i->s
                    - w(k,I_Praci_g) & ! [loss] i->g
                    - w(k,I_Pgaci  )   ! [loss] i->g
          qi_t(k) = max( qi_t(k), -w(k,I_dqi_dt) )
       end do

       do k = KS, KE
          qs_t(k) = + w(k,I_Psdep  ) & ! [prod] v->s
                    + w(k,I_Psacw  ) & ! [prod] c->s
                    + w(k,I_Psfw   ) & ! [prod] c->s
                    + w(k,I_Piacr_s) & ! [prod] r->s
                    + w(k,I_Psacr_s) & ! [prod] r->s
                    + w(k,I_Psaut  ) & ! [prod] i->s
                    + w(k,I_Praci_s) & ! [prod] i->s
                    + w(k,I_Psaci  ) & ! [prod] i->s
                    + w(k,I_Psfi   ) & ! [prod] i->s
                    - w(k,I_Pssub  ) & ! [loss] s->v
                    - w(k,I_Psmlt  ) & ! [loss] s->r
                    - w(k,I_Pgaut  ) & ! [loss] s->g
                    - w(k,I_Pracs  ) & ! [loss] s->g
                    - w(k,I_Pgacs  )   ! [loss] s->g
          qs_t(k) = max( qs_t(k), -w(k,I_dqs_dt) )
       end do

       do k = KS, KE
          qg_t(k) = + w(k,I_Pgdep  ) & ! [prod] v->g
                    + w(k,I_Pgacw  ) & ! [prod] c->g
                    + w(k,I_Piacr_g) & ! [prod] r->g
                    + w(k,I_Psacr_g) & ! [prod] r->g
                    + w(k,I_Pgacr  ) & ! [prod] r->g
                    + w(k,I_Pgfrz  ) & ! [prod] r->g
                    + w(k,I_Praci_g) & ! [prod] i->g
                    + w(k,I_Pgaci  ) & ! [prod] i->g
                    + w(k,I_Pgaut  ) & ! [prod] s->g
                    + w(k,I_Pracs  ) & ! [prod] s->g
                    + w(k,I_Pgacs  ) & ! [prod] s->g
                    - w(k,I_Pgsub  ) & ! [loss] g->v
                    - w(k,I_Pgmlt  )   ! [loss] g->r
          qg_t(k) = max( qg_t(k), -w(k,I_dqg_dt) )
       end do

       do k = KS, KE
          qv_t(k) = - ( qc_t(k) + qr_t(k) + qi_t(k) + qs_t(k) + qg_t(k) )
       end do

       do k = KS, KE
          QTRC0(k,i,j,I_QV) = QTRC0(k,i,j,I_QV) + qv_t(k) * dt
       end do
       do k = KS, KE
          QTRC0(k,i,j,I_QC) = QTRC0(k,i,j,I_QC) + qc_t(k) * dt
       end do
       do k = KS, KE
          QTRC0(k,i,j,I_QR) = QTRC0(k,i,j,I_QR) + qr_t(k) * dt
       end do
       do k = KS, KE
          QTRC0(k,i,j,I_QI) = QTRC0(k,i,j,I_QI) + qi_t(k) * dt
       end do
       do k = KS, KE
          QTRC0(k,i,j,I_QS) = QTRC0(k,i,j,I_QS) + qs_t(k) * dt
       end do
       do k = KS, KE
          QTRC0(k,i,j,I_QG) = QTRC0(k,i,j,I_QG) + qg_t(k) * dt
       end do
       do k = KS, KE
          cv_t = CV_VAPOR * qv_t(k) &
               + CV_WATER * ( qc_t(k) + qr_t(k) ) &
               + CV_ICE   * ( qi_t(k) + qs_t(k) + qg_t(k) )
          cvtot(k) = CVtot0(k,i,j) + cv_t * dt
       end do

       do k = KS, KE
          e_t = - LHV * qv_t(k) + LHF * ( qi_t(k) + qs_t(k) + qg_t(k) ) ! internal energy change
          RHOE_t(k,i,j) = dens(k) * e_t
          TEMP0(k,i,j) = ( temp(k) * CVtot0(k,i,j) + e_t * dt ) / cvtot(k)
       end do

       do k = KS, KE
          CVtot0(k,i,j) = cvtot(k)
       end do

       do k = KS, KE
          cp_t = CP_VAPOR * qv_t(k) &
               + CP_WATER * ( qc_t(k) + qr_t(k) ) &
               + CP_ICE   * ( qi_t(k) + qs_t(k) + qg_t(k) )
          cptot = CPtot0(k,i,j) + cp_t * dt
          CPtot0(k,i,j) = cptot
       end do

       if ( hist_flag ) then
          do ip = 1, w_nmax
             if ( HIST_sw(ip) ) w3d(k,i,j,ip) = w(k,ip)
          enddo
       end if

    enddo
    enddo

    do ip = 1, w_nmax
       if ( HIST_sw(ip) ) call FILE_HISTORY_put( HIST_id(ip), w3d(:,:,:,ip) )
    enddo

    call PROF_rapend  ('MP_tomita08', 3)

    return
  end subroutine MP_tomita08

  !-----------------------------------------------------------------------------
  !> Lin-type cold rain microphysics (terminal velocity)
!OCL SERIAL
  subroutine ATMOS_PHY_MP_tomita08_terminal_velocity( &
       KA, KS, KE, &
       DENS0, TEMP0, RHOQ0, &
       vterm                )
    use scale_const, only: &
       TEM00 => CONST_TEM00
    implicit none
    integer,  intent(in)  :: KA, KS, KE

    real(RP), intent(in)  :: DENS0(KA)
    real(RP), intent(in)  :: TEMP0(KA)
    real(RP), intent(in)  :: RHOQ0(KA,QA_MP-1)

    real(RP), intent(out) :: vterm(KA,QA_MP-1)

    real(RP) :: dens(KA)
    real(RP) :: temc(KA)
    real(RP) :: qr(KA), qi(KA), qs(KA), qg(KA)

    real(RP) :: rho_fact(KA) ! density factor

    real(RP) :: N0r(KA), N0s(KA), N0g(KA)
    real(RP) :: RLMDr, RLMDs, RLMDg
    real(RP) :: RLMDr_dr, RLMDs_ds, RLMDg_dg

    !---< Roh and Satoh (2014) >---
    real(RP) :: tems, Xs2
    real(RP) :: MOMs_0bs, RMOMs_Vt(KA)
    real(RP) :: coef_at(4), coef_bt(4)
    real(RP) :: loga_, b_, nm

    real(RP) :: zerosw
    integer  :: k
    !---------------------------------------------------------------------------

    do k = KS, KE
       ! store to work
       dens(k) = DENS0(k)
       temc(k) = TEMP0(k) - TEM00
       qr  (k) = RHOQ0(k,I_hyd_QR) / dens(k)
       qi  (k) = RHOQ0(k,I_hyd_QI) / dens(k)
       qs  (k) = RHOQ0(k,I_hyd_QS) / dens(k)
       qg  (k) = RHOQ0(k,I_hyd_QG) / dens(k)

       rho_fact(k) = sqrt( dens00 / dens(k) )
    end do

    if ( enable_WDXZ2014 ) then
       ! Wainwright et al. (2014)
       do k = KS, KE
          ! intercept parameter N0
          N0r(k) = 1.16E+5_RP * exp( log( max( dens(k)*qr(k)*1000.0_RP, 1.E-2_RP ) )*0.477_RP )
          N0s(k) = 4.58E+9_RP * exp( log( max( dens(k)*qs(k)*1000.0_RP, 1.E-2_RP ) )*0.788_RP )
          N0g(k) = 9.74E+8_RP * exp( log( max( dens(k)*qg(k)*1000.0_RP, 1.E-2_RP ) )*0.816_RP )
       end do
    else
       do k = KS, KE
          ! intercept parameter N0
          N0r(k) = N0r_def ! Marshall and Palmer (1948)
          N0s(k) = N0s_def ! Gunn and Marshall (1958)
          N0g(k) = N0g_def
       end do
    end if


    do k = KS, KE
       !---< terminal velocity >
       zerosw = 0.5_RP - sign(0.5_RP, qi(k) - 1.E-8_RP )
       vterm(k,I_hyd_QI) = -3.29_RP * exp( log( dens(k)*qi(k)+zerosw )*0.16_RP ) * ( 1.0_RP-zerosw )
    end do


    do k = KS, KE
       ! slope parameter lambda (Rain)
       zerosw = 0.5_RP - sign(0.5_RP, qr(k) - 1.E-12_RP )
       RLMDr  = sqrt(sqrt( dens(k) * qr(k) / ( Ar * N0r(k) * GAM_1br ) + zerosw )) * ( 1.0_RP-zerosw )
       RLMDr_dr = sqrt( RLMDr )       ! **Dr
       vterm(k,I_hyd_QR) = -Cr * rho_fact(k) * GAM_1brdr / GAM_1br * RLMDr_dr
    end do


    if ( enable_RS2014 ) then
       !---< modification by Roh and Satoh (2014) >---
       ! bimodal size distribution of snow
       do k = KS, KE
          zerosw = 0.5_RP - sign(0.5_RP, dens(k) * qs(k) - 1.E-12_RP )
          Xs2    = dens(k) * qs(k) / As

          tems       = min( -0.1_RP, temc(k) )
          coef_at(1) = coef_a01 + tems * ( coef_a02 + tems * ( coef_a05 + tems * coef_a09 ) )
          coef_at(2) = coef_a03 + tems * ( coef_a04 + tems *   coef_a07 )
          coef_at(3) = coef_a06 + tems *   coef_a08
          coef_at(4) = coef_a10
          coef_bt(1) = coef_b01 + tems * ( coef_b02 + tems * ( coef_b05 + tems * coef_b09 ) )
          coef_bt(2) = coef_b03 + tems * ( coef_b04 + tems *   coef_b07 )
          coef_bt(3) = coef_b06 + tems *   coef_b08
          coef_bt(4) = coef_b10
          ! 0 + Bs(=2) moment
          nm = 2.0_RP
          loga_  = coef_at(1) + nm * ( coef_at(2) + nm * ( coef_at(3) + nm * coef_at(4) ) )
          b_  = coef_bt(1) + nm * ( coef_bt(2) + nm * ( coef_bt(3) + nm * coef_bt(4) ) )
          MOMs_0bs = exp(ln10*loga_) * exp(log(Xs2+zerosw)*b_) * ( 1.0_RP-zerosw )
          ! Bs(=2) + Ds(=0.25) moment
          nm = 2.25_RP
          loga_  = coef_at(1) + nm * ( coef_at(2) + nm * ( coef_at(3) + nm * coef_at(4) ) )
          b_  = coef_bt(1) + nm * ( coef_bt(2) + nm * ( coef_bt(3) + nm * coef_bt(4) ) )
          RMOMs_Vt(k) = exp(ln10*loga_) * exp(log(Xs2+zerosw)*b_) * ( 1.0_RP-zerosw ) / ( MOMs_0bs + zerosw )
       end do
    else
       ! slope parameter lambda (Snow)
       do k = KS, KE
          zerosw = 0.5_RP - sign(0.5_RP, qs(k) - 1.E-12_RP )
          RLMDs  = sqrt(sqrt( dens(k) * qs(k) / ( As * N0s(k) * GAM_1bs ) + zerosw )) * ( 1.0_RP-zerosw )
          RLMDs_ds = sqrt( sqrt(RLMDs) ) ! **Ds
          RMOMs_Vt(k) = GAM_1bsds / GAM_1bs * RLMDs_ds
       end do
    end if

    do k = KS, KE
       vterm(k,I_hyd_QS) = -Cs * rho_fact(k) * RMOMs_Vt(k)
    end do


    do k = KS, KE
       ! slope parameter lambda (Graupel)
       zerosw = 0.5_RP - sign(0.5_RP, qg(k) - 1.E-12_RP )
       RLMDg  = sqrt(sqrt( dens(k) * qg(k) / ( Ag * N0g(k) * GAM_1bg ) + zerosw )) * ( 1.0_RP-zerosw )
       RLMDg_dg = sqrt( RLMDg )       ! **Dg
       vterm(k,I_hyd_QG) = -Cg * rho_fact(k) * GAM_1bgdg / GAM_1bg * RLMDg_dg
    enddo


!OCL XFILL
    do k = KS, KE
       vterm(k,I_hyd_QC) = 0.0_RP
    end do

    if ( nofall_qr ) then
!OCL XFILL
       do k = KS, KE
          vterm(k,I_hyd_QR) = 0.0_RP
       enddo
    endif

    if ( nofall_qi ) then
!OCL XFILL
       do k = KS, KE
          vterm(k,I_hyd_QI) = 0.0_RP
       enddo
    endif

    if ( nofall_qs ) then
!OCL XFILL
       do k = KS, KE
          vterm(k,I_hyd_QS) = 0.0_RP
       enddo
    endif

    if ( nofall_qg ) then
!OCL XFILL
       do k = KS, KE
          vterm(k,I_hyd_QG) = 0.0_RP
       enddo
    endif

    vterm(    1:KS-1,:) = 0.0_RP
    vterm(KE+1:KA   ,:) = 0.0_RP

    return
  end subroutine ATMOS_PHY_MP_tomita08_terminal_velocity

  !-----------------------------------------------------------------------------
  !> Calculate Cloud Fraction
  subroutine ATMOS_PHY_MP_tomita08_cloud_fraction( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       QTRC,           &
       mask_criterion, &
       cldfrac         )
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: QTRC(KA,IA,JA,QA_MP-1)
    real(RP), intent(in)  :: mask_criterion

    real(RP), intent(out) :: cldfrac(KA,IA,JA)

    real(RP) :: qhydro
    integer  :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do OMP_SCHEDULE_ &
    !$omp private(qhydro)
    do j  = JS, JE
    do i  = IS, IE
    do k  = KS, KE
       qhydro = QTRC(k,i,j,I_hyd_QC) + QTRC(k,i,j,I_hyd_QR) &
              + QTRC(k,i,j,I_hyd_QI) + QTRC(k,i,j,I_hyd_QS) + QTRC(k,i,j,I_hyd_QG)
       cldfrac(k,i,j) = 0.5_RP + sign(0.5_RP,qhydro-mask_criterion)
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_MP_tomita08_cloud_fraction

  !-----------------------------------------------------------------------------
  !> Calculate Effective Radius
  subroutine ATMOS_PHY_MP_tomita08_effective_radius( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       DENS0, TEMP0, QTRC0, &
       Re                   )
    use scale_const, only: &
       PI    => CONST_PI,     &
       dens_w => CONST_DWATR, &
       dens_i => CONST_DICE,  &
       TEM00 => CONST_TEM00
    use scale_atmos_hydrometeor, only: &
       N_HYD, &
       I_HC, &
       I_HR, &
       I_HI, &
       I_HS, &
       I_HG
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: DENS0(KA,IA,JA) !> density                   [kg/m3]
    real(RP), intent(in)  :: TEMP0(KA,IA,JA) !> temperature               [K]
    real(RP), intent(in)  :: QTRC0(KA,IA,JA,QA_MP-1)

    real(RP), intent(out) :: Re   (KA,IA,JA,N_HYD) ! effective radius          [cm]
    real(RP) :: dens(KA)
    real(RP) :: temc(KA)
    real(RP) :: qr(KA), qs(KA), qg(KA)
    real(RP) :: Nc(KA)              !< Number concentration of cloud water [1/m3]
    real(RP) :: N0r(KA), N0s(KA), N0g(KA)
    real(RP) :: RLMDr, RLMDs, RLMDg

    real(RP), parameter :: um2cm = 100.0_RP

    !---< Roh and Satoh (2014) >---
    real(RP) :: tems, Xs2
    real(RP) :: coef_at(4), coef_bt(4)
    real(RP) :: loga_, b_, nm

    real(RP) :: zerosw
    integer  :: k, i, j, ih
    !---------------------------------------------------------------------------

    !$omp parallel do OMP_SCHEDULE_
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       Re(k,i,j,I_HI) =  re_qi * um2cm
    end do
    end do
    end do
    !$omp parallel do OMP_SCHEDULE_
    do ih = I_HG+1, N_HYD
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       Re(k,i,j,ih) = 0.0_RP
    end do
    end do
    end do
    end do

    if ( const_rec .or. fixed_re ) then

       !$omp parallel do OMP_SCHEDULE_
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          Re(k,i,j,I_HC) = re_qc * um2cm
       end do
       end do
       end do

    else

       !$omp parallel do OMP_SCHEDULE_ &
       !$omp private(Nc)
       do j = JS, JE
       do i = IS, IE
          if ( do_couple_aerosol ) then
             do k = KS, KE
                ! Nc(k) = max( CCN(k,i,j), Nc_def(i,j)*1.E+6_RP ) ! [#/m3] tentatively off the CCN effect
                Nc(k) = Nc_def(i,j) * 1.E+6_RP ! [#/m3]
             enddo
          else
             do k = KS, KE
                Nc(k) = Nc_def(i,j) * 1.E+6_RP ! [#/m3]
             enddo
          endif

          do k  = KS, KE
             Re(k,i,j,I_HC) = 1.1_RP &
                            * ( DENS0(k,i,j) * QTRC0(k,i,j,I_hyd_QC) / Nc(k) / ( 4.0_RP / 3.0_RP * PI * dens_w ) )**(1.0_RP/3.0_RP)
             Re(k,i,j,I_HC) = min( 1.E-3_RP, max( 1.E-6_RP, Re(k,i,j,I_HC) ) ) * um2cm
          enddo
       enddo
       enddo

    endif

    if ( fixed_re ) then

       !$omp parallel do OMP_SCHEDULE_
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          Re(k,i,j,I_HR) =  10000.E-6_RP * um2cm
       end do
       end do
       end do
       !$omp parallel do OMP_SCHEDULE_
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          Re(k,i,j,I_HS) =  re_qi * um2cm
       end do
       end do
       end do
       !$omp parallel do OMP_SCHEDULE_
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
       Re(k,i,j,I_HG) =  re_qi * um2cm
       end do
       end do
       end do

    else

#ifndef __GFORTRAN__
       !$omp parallel do default(none)                                                          &
       !$omp shared(JS,JE,IS,IE,KS,KE,DENS0,TEMP0,QTRC0,enable_WDXZ2014,N0r_def)                &
       !$omp shared(N0s_def,N0g_def,Ar,GAM_1br,Re,As,GAM_1bs)                                   &
       !$omp shared(enable_RS2014,ln10,Ag,GAM_1bg)                                              &
       !$omp private(i,j,k,dens,temc,qr,qs,qg,N0r,N0s,N0g,zerosw,RLMDr,RLMDs,Xs2,nm,tems,loga_) &
       !$omp private(b_,RLMDg,coef_at,coef_bt) &
       !$omp OMP_SCHEDULE_ collapse(2)
#else
       !$omp parallel do default(shared)                                                        &
       !$omp private(i,j,k,dens,temc,qr,qs,qg,N0r,N0s,N0g,zerosw,RLMDr,RLMDs,Xs2,nm,tems,loga_) &
       !$omp private(b_,RLMDg,coef_at,coef_bt) &
       !$omp OMP_SCHEDULE_ collapse(2)
#endif
       do j  = JS, JE
       do i  = IS, IE

          do k = KS, KE
             dens(k) = DENS0(k,i,j)
             temc(k) = TEMP0(k,i,j) - TEM00
             qr  (k) = QTRC0(k,i,j,I_hyd_QR)
             qs  (k) = QTRC0(k,i,j,I_hyd_QS)
             qg  (k) = QTRC0(k,i,j,I_hyd_QG)
          end do

          ! intercept parameter N0
          if ( enable_WDXZ2014 ) then
             ! Wainwright et al. (2014)
             do k = KS, KE
                N0r(k) = 1.16E+5_RP * exp( log( max( dens(k)*qr(k)*1000.0_RP, 1.E-2_RP ) )*0.477_RP )
                N0s(k) = 4.58E+9_RP * exp( log( max( dens(k)*qs(k)*1000.0_RP, 1.E-2_RP ) )*0.788_RP )
                N0g(k) = 9.74E+8_RP * exp( log( max( dens(k)*qg(k)*1000.0_RP, 1.E-2_RP ) )*0.816_RP )
             end do
          else
             do k = KS, KE
                N0r(k) = N0r_def ! Marshall and Palmer (1948)
                N0s(k) = N0s_def ! Gunn and Marshall (1958)
                N0g(k) = N0g_def
             end do
          end if

          ! slope parameter lambda
          do k = KS, KE
             zerosw = 0.5_RP - sign(0.5_RP, qr(k) - 1.E-12_RP )
             RLMDr = sqrt(sqrt( dens(k) * qr(k) / ( Ar * N0r(k) * GAM_1br ) + zerosw )) * ( 1.0_RP-zerosw )
             ! Effective radius is defined by r3m/r2m = 1.5/lambda
             Re(k,i,j,I_HR) = 1.5_RP * RLMDr * um2cm
          end do

          if ( enable_RS2014 ) then
             do k = KS, KE
                zerosw = 0.5_RP - sign(0.5_RP, qs(k) - 1.E-12_RP )
                !---< modification by Roh and Satoh (2014) >---
                ! bimodal size distribution of snow
                zerosw = 0.5_RP - sign(0.5_RP, dens(k) * qs(k) - 1.E-12_RP )
                Xs2    = dens(k) * qs(k) / As

                tems       = min( -0.1_RP, temc(k) )
                coef_at(1) = coef_a01 + tems * ( coef_a02 + tems * ( coef_a05 + tems * coef_a09 ) )
                coef_at(2) = coef_a03 + tems * ( coef_a04 + tems *   coef_a07 )
                coef_at(3) = coef_a06 + tems *   coef_a08
                coef_at(4) = coef_a10
                coef_bt(1) = coef_b01 + tems * ( coef_b02 + tems * ( coef_b05 + tems * coef_b09 ) )
                coef_bt(2) = coef_b03 + tems * ( coef_b04 + tems *   coef_b07 )
                coef_bt(3) = coef_b06 + tems *   coef_b08
                coef_bt(4) = coef_b10
                ! 1 + Bs(=2) moment
                nm = 3.0_RP
                loga_  = coef_at(1) + nm * ( coef_at(2) + nm * ( coef_at(3) + nm * coef_at(4) ) )
                   b_  = coef_bt(1) + nm * ( coef_bt(2) + nm * ( coef_bt(3) + nm * coef_bt(4) ) )

                Re(k,i,j,I_HS) = 0.5_RP * exp(ln10*loga_) * exp(log(Xs2+zerosw)*b_) * ( 1.0_RP-zerosw ) / ( Xs2+zerosw ) * um2cm
             end do
          else
             do k = KS, KE
                zerosw = 0.5_RP - sign(0.5_RP, qs(k) - 1.E-12_RP )
                RLMDs = sqrt(sqrt( dens(k) * qs(k) / ( As * N0s(k) * GAM_1bs ) + zerosw )) * ( 1.0_RP-zerosw )
                Re(k,i,j,I_HS) = 1.5_RP * RLMDs * um2cm
                ! Re(k,i,j,I_HS) = dens * qs / N0s / ( 2.0_RP / 3.0_RP * PI * dens_i ) / RLMDs**3 * um2cm
             end do
          end if

          do k = KS, KE
             zerosw = 0.5_RP - sign(0.5_RP, qg(k) - 1.E-12_RP )
             RLMDg = sqrt(sqrt( dens(k) * qg(k) / ( Ag * N0g(k) * GAM_1bg ) + zerosw )) * ( 1.0_RP-zerosw )
             Re(k,i,j,I_HG) = 1.5_RP * RLMDg * um2cm
             ! Re(k,i,j,I_HG) = dens * qg / N0g / ( 2.0_RP / 3.0_RP * PI * dens_i ) / RLMDg**3 * um2cm
          enddo

       enddo
       enddo

    endif

    return
  end subroutine ATMOS_PHY_MP_tomita08_effective_radius

  !-----------------------------------------------------------------------------
!OCL SERIAL
  subroutine MP_tomita08_BergeronParam( &
       KA, KS, KE, &
       temp,       &
       a1, a2, ma2 )
    use scale_const, only: &
       TEM00 => CONST_TEM00
    implicit none
    integer, intent(in) :: KA, KS, KE

    real(RP), intent(in)  :: temp(KA)
    real(RP), intent(out) :: a1(KA)
    real(RP), intent(out) :: a2(KA)
    real(RP), intent(out) :: ma2(KA)

    real(RP), parameter :: a1_tab(32) = (/ &
         0.0001E-7_RP, 0.7939E-7_RP, 0.7841E-6_RP, 0.3369E-5_RP, 0.4336E-5_RP, &
         0.5285E-5_RP, 0.3728E-5_RP, 0.1852E-5_RP, 0.2991E-6_RP, 0.4248E-6_RP, &
         0.7434E-6_RP, 0.1812E-5_RP, 0.4394E-5_RP, 0.9145E-5_RP, 0.1725E-4_RP, &
         0.3348E-4_RP, 0.1725E-4_RP, 0.9175E-5_RP, 0.4412E-5_RP, 0.2252E-5_RP, &
         0.9115E-6_RP, 0.4876E-6_RP, 0.3473E-6_RP, 0.4758E-6_RP, 0.6306E-6_RP, &
         0.8573E-6_RP, 0.7868E-6_RP, 0.7192E-6_RP, 0.6513E-6_RP, 0.5956E-6_RP, &
         0.5333E-6_RP, 0.4834E-6_RP  /)
    real(RP), parameter :: a2_tab(32) = (/ &
         0.0100_RP, 0.4006_RP, 0.4831_RP, 0.5320_RP, 0.5307_RP, &
         0.5319_RP, 0.5249_RP, 0.4888_RP, 0.3849_RP, 0.4047_RP, &
         0.4318_RP, 0.4771_RP, 0.5183_RP, 0.5463_RP, 0.5651_RP, &
         0.5813_RP, 0.5655_RP, 0.5478_RP, 0.5203_RP, 0.4906_RP, &
         0.4447_RP, 0.4126_RP, 0.3960_RP, 0.4149_RP, 0.4320_RP, &
         0.4506_RP, 0.4483_RP, 0.4460_RP, 0.4433_RP, 0.4413_RP, &
         0.4382_RP, 0.4361_RP  /)

    real(RP) :: temc
    integer  :: itemc
    real(RP) :: fact

    integer :: k
    !---------------------------------------------------------------------------

    do k = KS, KE
       temc  = min( max( temp(k)-TEM00, -30.99_RP ), 0.0_RP ) ! 0C <= T  <  31C
       itemc = int( -temc ) + 1                            ! 1  <= iT <= 31
       fact  = - ( temc + real(itemc-1,kind=8) )

       a1(k) = ( 1.0_RP-fact ) * a1_tab(itemc  ) &
             + (        fact ) * a1_tab(itemc+1)

       a2(k) = ( 1.0_RP-fact ) * a2_tab(itemc  ) &
             + (        fact ) * a2_tab(itemc+1)

       ma2(k) = 1.0_RP - a2(k)

       a1(k) = a1(k) * exp( log(1.E-3_RP)*ma2(k) ) ! [g->kg]

    end do

    return
  end subroutine MP_tomita08_BergeronParam

  !-----------------------------------------------------------------------------
  !> Calculate mass ratio of each category
  subroutine ATMOS_PHY_MP_tomita08_qtrc2qhyd( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       QTRC, &
       Qe    )
    use scale_atmos_hydrometeor, only: &
       N_HYD, &
       I_HC, &
       I_HR, &
       I_HI, &
       I_HS, &
       I_HG
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: QTRC(KA,IA,JA,QA_MP-1)

    real(RP), intent(out) :: Qe(KA,IA,JA,N_HYD) ! mass ratio of each cateory [kg/kg]

    integer :: k, i, j, ih
    !---------------------------------------------------------------------------

    !$omp parallel do OMP_SCHEDULE_
!OCL XFILL
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       Qe(k,i,j,I_HC) = QTRC(k,i,j,I_hyd_QC)
    end do
    end do
    end do
    !$omp parallel do OMP_SCHEDULE_
!OCL XFILL
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       Qe(k,i,j,I_HR) = QTRC(k,i,j,I_hyd_QR)
    end do
    end do
    end do
    !$omp parallel do OMP_SCHEDULE_
!OCL XFILL
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       Qe(k,i,j,I_HI) = QTRC(k,i,j,I_hyd_QI)
    end do
    end do
    end do
    !$omp parallel do OMP_SCHEDULE_
!OCL XFILL
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       Qe(k,i,j,I_HS) = QTRC(k,i,j,I_hyd_QS)
    end do
    end do
    end do
    !$omp parallel do OMP_SCHEDULE_
!OCL XFILL
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       Qe(k,i,j,I_HG) = QTRC(k,i,j,I_hyd_QG)
    end do
    end do
    end do
    !$omp parallel do OMP_SCHEDULE_
!OCL XFILL
    do ih = I_HG+1, N_HYD
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       Qe(k,i,j,ih) = 0.0_RP
    end do
    end do
    end do
    end do

    return
  end subroutine ATMOS_PHY_MP_tomita08_qtrc2qhyd

  !-----------------------------------------------------------------------------
  !> get mass ratio of each category
  subroutine ATMOS_PHY_MP_tomita08_qhyd2qtrc( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       Qe,  &
       QTRC )
    use scale_atmos_hydrometeor, only: &
       N_HYD, &
       I_HC,  &
       I_HR,  &
       I_HI,  &
       I_HS,  &
       I_HG,  &
       I_HH
    implicit none

    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE
    real(RP), intent(in)  :: Qe  (KA,IA,JA,N_HYD)   ! mass ratio of each cateory [kg/kg]
    real(RP), intent(out) :: QTRC(KA,IA,JA,QA_MP-1)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do OMP_SCHEDULE_
!OCL XFILL
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       QTRC(k,i,j,I_hyd_QC) = Qe(k,i,j,I_HC)
    end do
    end do
    end do
    !$omp parallel do OMP_SCHEDULE_
!OCL XFILL
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       QTRC(k,i,j,I_hyd_QR) = Qe(k,i,j,I_HR)
    end do
    end do
    end do
    !$omp parallel do OMP_SCHEDULE_
!OCL XFILL
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       QTRC(k,i,j,I_hyd_QI) = Qe(k,i,j,I_HI)
    end do
    end do
    end do
    !$omp parallel do OMP_SCHEDULE_
!OCL XFILL
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       QTRC(k,i,j,I_hyd_QS) = Qe(k,i,j,I_HS)
    end do
    end do
    end do
    !$omp parallel do OMP_SCHEDULE_
!OCL XFILL
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       QTRC(k,i,j,I_hyd_QG) = Qe(k,i,j,I_HG) + Qe(k,i,j,I_HH)
    end do
    end do
    end do

    return
  end subroutine ATMOS_PHY_MP_tomita08_qhyd2qtrc

end module scale_atmos_phy_mp_tomita08
