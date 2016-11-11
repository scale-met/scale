!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics
!!
!! @par Description
!!          Cloud Microphysics by Lin-type parametarization
!!          Reference: Tomita(2008)
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-03-25 (H.Yashiro)  [new]
!! @li      2015-09-08 (Y.Sato)     [add] Add evaporated cloud number concentration
!!
!<
!-------------------------------------------------------------------------------
module scale_atmos_phy_mp_tomita08
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index

  use scale_atmos_hydrometeor, only: &
     N_HYD, &
     I_QV, &
     I_QC, &
     I_QR, &
     I_QI, &
     I_QS, &
     I_QG, &
     I_HC, &
     I_HR, &
     I_HI, &
     I_HS, &
     I_HG
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_MP_tomita08_config
  public :: ATMOS_PHY_MP_tomita08_setup
  public :: ATMOS_PHY_MP_tomita08
  public :: ATMOS_PHY_MP_tomita08_CloudFraction
  public :: ATMOS_PHY_MP_tomita08_EffectiveRadius
  public :: ATMOS_PHY_MP_tomita08_Mixingratio

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, private, parameter :: QA_MP  = 6

  character(len=H_SHORT), public, target :: ATMOS_PHY_MP_tomita08_NAME(QA_MP)
  character(len=H_MID)  , public, target :: ATMOS_PHY_MP_tomita08_DESC(QA_MP)
  character(len=H_SHORT), public, target :: ATMOS_PHY_MP_tomita08_UNIT(QA_MP)

  real(RP), public, target :: ATMOS_PHY_MP_tomita08_DENS(N_HYD) ! hydrometeor density [kg/m3]=[g/L]

  data ATMOS_PHY_MP_tomita08_NAME / &
                 'QV', &
                 'QC', &
                 'QR', &
                 'QI', &
                 'QS', &
                 'QG'  /

  data ATMOS_PHY_MP_tomita08_DESC / &
                 'Ratio of Water Vapor mass to total mass (Specific humidity)',   &
                 'Ratio of Cloud Water mass to total mass',   &
                 'Ratio of Rain Water mass to total mass',    &
                 'Ratio of Cloud Ice mixing ratio to total mass',     &
                 'Ratio of Snow mixing ratio to total mass',          &
                 'Ratio of Graupel mixing ratio to total mass'        /

  data ATMOS_PHY_MP_tomita08_UNIT / &
                 'kg/kg',  &
                 'kg/kg',  &
                 'kg/kg',  &
                 'kg/kg',  &
                 'kg/kg',  &
                 'kg/kg'   /

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: MP_tomita08
  private :: MP_tomita08_vterm
  private :: MP_tomita08_BergeronParam

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,  private, parameter :: I_mp_QC =  1
  integer,  private, parameter :: I_mp_QR =  2
  integer,  private, parameter :: I_mp_QI =  3
  integer,  private, parameter :: I_mp_QS =  4
  integer,  private, parameter :: I_mp_QG =  5
  integer,  private            :: QS_MP
  integer,  private            :: QE_MP

  logical,  private            :: MP_donegative_fixer  = .true.  ! apply negative fixer?
  logical,  private            :: MP_doprecipitation   = .true.  ! apply sedimentation (precipitation)?
  logical,  private            :: MP_couple_aerosol    = .false. ! apply CCN effect?
  real(RP), private            :: MP_limit_negative    = 1.0_RP  ! Abort if abs(fixed negative vaue) > abs(MP_limit_negative)

  logical,  private            :: MP_doexpricit_icegen = .false. ! apply explicit ice generation?
  logical,  private            :: only_liquid          = .false.
  real(RP), private            :: sw_useicegen         = 0.0_RP

  real(RP), private            :: dens00 = 1.28_RP !< standard density [kg/m3]

  ! Parameter for Marshall-Palmer distribution
  real(RP), private            :: N0r = 8.E6_RP !< intercept parameter for rain    [1/m4]
  real(RP), private            :: N0s = 3.E6_RP !< intercept parameter for snow    [1/m4]
  real(RP), private            :: N0g = 4.E6_RP !< intercept parameter for graupel [1/m4]

  real(RP), private            :: dens_s = 100.0_RP !< density of snow    [kg/m3]
  real(RP), private            :: dens_g = 400.0_RP !< density of graupel [kg/m3]
                                                    !   graupel : 400
                                                    !   hail    : 917
  real(RP), private            :: drag_g = 0.6_RP   !< drag coefficient for graupel
  real(RP), private            :: re_qc  =  8.E-6_RP ! effective radius for cloud water
  real(RP), private            :: re_qi  = 40.E-6_RP ! effective radius for cloud ice
  real(RP), private            :: Cr     = 130.0_RP
  real(RP), private            :: Cs     = 4.84_RP

  ! Empirical parameter
  real(RP), private            :: Ar, As, Ag
  real(RP), private            :: Br, Bs, Bg
  real(RP), private            :: Cg
  real(RP), private            :: Dr, Ds, Dg

  ! GAMMA function
  real(RP), private            :: GAM, GAM_2, GAM_3

  real(RP), private            :: GAM_1br, GAM_2br, GAM_3br
  real(RP), private            :: GAM_3dr
  real(RP), private            :: GAM_6dr
  real(RP), private            :: GAM_1brdr
  real(RP), private            :: GAM_5dr_h

  real(RP), private            :: GAM_1bs, GAM_2bs, GAM_3bs
  real(RP), private            :: GAM_3ds
  real(RP), private            :: GAM_1bsds
  real(RP), private            :: GAM_5ds_h

  real(RP), private            :: GAM_1bg, GAM_3dg
  real(RP), private            :: GAM_1bgdg
  real(RP), private            :: GAM_5dg_h

  !---< Khairoutdinov and Kogan (2000) >---
  real(RP), private            :: sw_kk2000  = 0.0_RP !< switch for k-k scheme

  !---< Roh and Satoh (2014) >---
  real(RP), private            :: sw_roh2014 = 0.0_RP !< switch for Roh scheme
  real(RP), private, parameter :: coef_a(10) = (/ 5.065339_RP, -0.062659_RP, -3.032362_RP, 0.029469_RP, -0.000285_RP, &
                                                  0.31255_RP,   0.000204_RP,  0.003199_RP, 0.0_RP,      -0.015952_RP  /)
  real(RP), private, parameter :: coef_b(10) = (/ 0.476221_RP, -0.015896_RP,  0.165977_RP, 0.007468_RP, -0.000141_RP, &
                                                  0.060366_RP,  0.000079_RP,  0.000594_RP, 0.0_RP,      -0.003577_RP  /)

  ! Accretion parameter
  real(RP), private            :: Eiw = 1.0_RP  !< collection efficiency of cloud ice for cloud water
  real(RP), private            :: Erw = 1.0_RP  !< collection efficiency of rain    for cloud water
  real(RP), private            :: Esw = 1.0_RP  !< collection efficiency of snow    for cloud water
  real(RP), private            :: Egw = 1.0_RP  !< collection efficiency of graupel for cloud water
  real(RP), private            :: Eri = 1.0_RP  !< collection efficiency of rain    for cloud ice
  real(RP), private            :: Esi = 1.0_RP  !< collection efficiency of snow    for cloud ice
  real(RP), private            :: Egi = 0.1_RP  !< collection efficiency of graupel for cloud ice
  real(RP), private            :: Esr = 1.0_RP  !< collection efficiency of snow    for rain
  real(RP), private            :: Egr = 1.0_RP  !< collection efficiency of graupel for rain
  real(RP), private            :: Egs = 1.0_RP  !< collection efficiency of graupel for snow
  real(RP), private            :: gamma_sacr = 25.E-3_RP !< effect of low temperature for Esi
  real(RP), private            :: gamma_gacs = 90.E-3_RP !< effect of low temperature for Egs

  ! Auto-conversion parameter
  real(RP), private, parameter   :: Nc_lnd = 2000.0_RP !< number concentration of cloud water (land)  [1/cc]
  real(RP), private, parameter   :: Nc_ocn =   50.0_RP !< number concentration of cloud water (ocean) [1/cc]
  real(RP), private, allocatable :: Nc_def(:,:)        !< number concentration of cloud water         [1/cc]

  real(RP), private            :: beta_saut  =  6.E-3_RP  !< auto-conversion factor beta  for ice
  real(RP), private            :: gamma_saut = 60.E-3_RP  !< auto-conversion factor gamma for ice
  real(RP), private            :: beta_gaut  =  1.E-3_RP  !< auto-conversion factor beta  for snow
  real(RP), private            :: gamma_gaut = 90.E-3_RP  !< auto-conversion factor gamma for snow
  real(RP), private            :: qicrt_saut =  0.0_RP    !< mixing ratio threshold for Psaut [kg/kg]
  real(RP), private            :: qscrt_gaut =  6.E-4_RP  !< mixing ratio threshold for Pgaut [kg/kg]

  ! Evaporation, Sublimation parameter
  real(RP), private, parameter :: Da0    = 2.428E-2_RP !< thermal diffusion coefficient of air at 0C,1atm [J/m/s/K]
  real(RP), private, parameter :: dDa_dT =  7.47E-5_RP !< Coefficient of Da depending on temperature      [J/m/s/K/K]
  real(RP), private, parameter :: Dw0    = 2.222E-5_RP !< diffusion coefficient of water vapor in the air at 0C,1atm [m2/s]
  real(RP), private, parameter :: dDw_dT =  1.37E-7_RP !< Coefficient of Dw depending on temperature                 [m2/s/K]
  real(RP), private, parameter :: mu0    = 1.718E-5_RP !< kinematic viscosity of air at 0C,1atm      [m2/s*kg/m3]
  real(RP), private, parameter :: dmu_dT =  5.28E-8_RP !< Coefficient of mu depending on temperature [m2/s/K*kg/m3]

  real(RP), private            :: f1r = 0.78_RP  !< ventilation factor 1 for rain
  real(RP), private            :: f2r = 0.27_RP  !< ventilation factor 2 for rain
  real(RP), private            :: f1s = 0.65_RP  !< ventilation factor 1 for snow
  real(RP), private            :: f2s = 0.39_RP  !< ventilation factor 2 for snow
  real(RP), private            :: f1g = 0.78_RP  !< ventilation factor 1 for graupel
  real(RP), private            :: f2g = 0.27_RP  !< ventilation factor 2 for graupel

  ! Freezing parameter
  real(RP), private            :: A_frz = 0.66_RP  !< freezing factor [/K]
  real(RP), private            :: B_frz = 100.0_RP !< freezing factor [/m3/s]

  ! Bergeron process parameter
  real(RP), private            :: mi40  = 2.46E-10_RP !< mass              of a 40 micron ice crystal [kg]
  real(RP), private            :: mi50  = 4.80E-10_RP !< mass              of a 50 micron ice crystal [kg]
  real(RP), private            :: vti50 = 1.0_RP      !< terminal velocity of a 50 micron ice crystal [m/s]
  real(RP), private            :: Ri50  = 5.E-5_RP    !< radius            of a 50 micron ice crystal [m]

  integer,  private, parameter :: w_nmax = 42
  integer,  private, parameter :: I_dqv_dt  =  1 !
  integer,  private, parameter :: I_dqc_dt  =  2 !
  integer,  private, parameter :: I_dqr_dt  =  3 !
  integer,  private, parameter :: I_dqi_dt  =  4 !
  integer,  private, parameter :: I_dqs_dt  =  5 !
  integer,  private, parameter :: I_dqg_dt  =  6 !
  integer,  private, parameter :: I_delta1  =  7 ! separation switch for r->s,g
  integer,  private, parameter :: I_delta2  =  8 ! separation switch for s->g
  integer,  private, parameter :: I_delta3  =  9 ! separation switch for ice sublimation
  integer,  private, parameter :: I_RLMDr   = 10
  integer,  private, parameter :: I_RLMDs   = 11
  integer,  private, parameter :: I_RLMDg   = 12
  integer,  private, parameter :: I_Piacr   = 13 ! r->s,g
  integer,  private, parameter :: I_Psacr   = 14 ! r->s,g
  integer,  private, parameter :: I_Praci   = 15 ! i->s,g
  integer,  private, parameter :: I_Psmlt   = 16 ! s->r
  integer,  private, parameter :: I_Pgmlt   = 17 ! g->r
  integer,  private, parameter :: I_Praut   = 18 ! c->r
  integer,  private, parameter :: I_Pracw   = 19 ! c->r
  integer,  private, parameter :: I_Psacw   = 20 ! c->s
  integer,  private, parameter :: I_Psfw    = 21 ! c->s
  integer,  private, parameter :: I_Pgacw   = 22 ! c->g
  integer,  private, parameter :: I_Prevp   = 23 ! r->v
  integer,  private, parameter :: I_Piacr_s = 24 ! r->s
  integer,  private, parameter :: I_Psacr_s = 25 ! r->s
  integer,  private, parameter :: I_Piacr_g = 26 ! r->g
  integer,  private, parameter :: I_Psacr_g = 27 ! r->g
  integer,  private, parameter :: I_Pgacr   = 28 ! r->g
  integer,  private, parameter :: I_Pgfrz   = 29 ! r->g
  integer,  private, parameter :: I_Psaut   = 30 ! i->s
  integer,  private, parameter :: I_Praci_s = 31 ! i->s
  integer,  private, parameter :: I_Psaci   = 32 ! i->s
  integer,  private, parameter :: I_Psfi    = 33 ! i->s
  integer,  private, parameter :: I_Praci_g = 34 ! i->g
  integer,  private, parameter :: I_Pgaci   = 35 ! i->g
  integer,  private, parameter :: I_Psdep   = 36 ! v->s
  integer,  private, parameter :: I_Pssub   = 37 ! s->v
  integer,  private, parameter :: I_Pgaut   = 38 ! s->g
  integer,  private, parameter :: I_Pracs   = 39 ! s->g
  integer,  private, parameter :: I_Pgacs   = 40 ! s->g
  integer,  private, parameter :: I_Pgdep   = 41 ! v->g
  integer,  private, parameter :: I_Pgsub   = 42 ! g->v

  character(len=H_SHORT), private :: w_name(w_nmax)

  data w_name / 'dqv_dt ', &
                'dqc_dt ', &
                'dqr_dt ', &
                'dqi_dt ', &
                'dqs_dt ', &
                'dqg_dt ', &
                'delta1 ', &
                'delta2 ', &
                'delta3 ', &
                'RLMDr  ', &
                'RLMDs  ', &
                'RLMDg  ', &
                'Piacr  ', &
                'Psacr  ', &
                'Praci  ', &
                'Psmlt  ', &
                'Pgmlt  ', &
                'Praut  ', &
                'Pracw  ', &
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
                'Psaut  ', &
                'Praci_s', &
                'Psaci  ', &
                'Psfi   ', &
                'Praci_g', &
                'Pgaci  ', &
                'Psdep  ', &
                'Pssub  ', &
                'Pgaut  ', &
                'Pracs  ', &
                'Pgacs  ', &
                'Pgdep  ', &
                'Pgsub  '  /

  real(RP), private, allocatable :: w3d(:,:,:,:) !< for history output

  integer,  private :: MP_ntmax_sedimentation = 1 ! number of time step for sedimentation

  integer,  private :: MP_NSTEP_SEDIMENTATION
  real(RP), private :: MP_RNSTEP_SEDIMENTATION
  real(DP), private :: MP_DTSEC_SEDIMENTATION

  logical, private :: debug

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Config
  subroutine ATMOS_PHY_MP_tomita08_config( &
       MP_TYPE, &
       QA, QS   )
    use scale_process, only: &
       PRC_MPIstop
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_regist
    implicit none

    character(len=*), intent(in) :: MP_TYPE
    integer, intent(out) :: QA
    integer, intent(out) :: QS
    !---------------------------------------------------------------------------

    if ( MP_TYPE /= 'TOMITA08' ) then
       write(*,*) 'xxx ATMOS_PHY_MP_TYPE is not TOMITA08. Check!'
       call PRC_MPIstop
    endif

    call ATMOS_HYDROMETEOR_regist( QS,                         & ! (out)
                                   1, 2, 3,                    & ! (in)
                                   ATMOS_PHY_MP_tomita08_NAME, & ! (in)
                                   ATMOS_PHY_MP_tomita08_DESC, & ! (in)
                                   ATMOS_PHY_MP_tomita08_UNIT  ) ! (in)

    QA = QA_MP
    QS_MP = QS
    QE_MP = QS + QA - 1

    I_QV = QS
    I_QC = QS + I_mp_QC
    I_QR = QS + I_mp_QR
    I_QI = QS + I_mp_QI
    I_QS = QS + I_mp_QS
    I_QG = QS + I_mp_QG

    return
  end subroutine ATMOS_PHY_MP_tomita08_config

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_MP_tomita08_setup
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       UNDEF  => CONST_UNDEF, &
       PI     => CONST_PI,    &
       GRAV   => CONST_GRAV,  &
       dens_w => CONST_DWATR,  &
       dens_i => CONST_DICE
    use scale_specfunc, only: &
       SF_gamma
    use scale_time, only: &
       TIME_DTSEC_ATMOS_PHY_MP
    use scale_grid, only: &
       CDZ => GRID_CDZ
    implicit none

    real(RP) :: autoconv_nc    = Nc_ocn  !< number concentration of cloud water [1/cc]
    logical  :: enable_kk2000  = .false. !< use scheme by Khairoutdinov and Kogan (2000)
    logical  :: enable_roh2014 = .false. !< use scheme by Roh and Satoh (2014)

    NAMELIST / PARAM_ATMOS_PHY_MP / &
       MP_doprecipitation,     &
       MP_donegative_fixer,    &
       MP_limit_negative,      &
       MP_doexpricit_icegen,   &
       MP_ntmax_sedimentation, &
       MP_couple_aerosol

    NAMELIST / PARAM_ATMOS_PHY_MP_TOMITA08 / &
       autoconv_nc,    &
       enable_kk2000,  &
       enable_roh2014, &
       dens_s,         &
       dens_g,         &
       re_qc,          &
       re_qi,          &
       Cr,             &
       Cs,             &
       drag_g,         &
       beta_saut,      &
       gamma_saut,     &
       gamma_sacr,     &
       N0s,            &
       N0g,            &
       debug

    real(RP), parameter :: max_term_vel = 10.0_RP  !-- terminal velocity for calculate dt of sedimentation
    integer :: nstep_max

    integer :: ierr
    integer :: i, j, ip
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[Cloud Microphysics] / Categ[ATMOS PHYSICS] / Origin[SCALElib]'
    if( IO_L ) write(IO_FID_LOG,*) '*** TOMITA08: 1-moment bulk 6 category'

    allocate( w3d(KA,IA,JA,w_nmax) )
    w3d(:,:,:,:) = 0.0_RP

    allocate( Nc_def(IA,JA) )

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_MP,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_MP. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_MP)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Enable negative fixer?                    : ', MP_donegative_fixer
    if( IO_L ) write(IO_FID_LOG,*) '*** Value limit of negative fixer (abs)       : ', abs(MP_limit_negative)
    if( IO_L ) write(IO_FID_LOG,*) '*** Enable sedimentation (precipitation)?     : ', MP_doprecipitation
    if( IO_L ) write(IO_FID_LOG,*) '*** Enable explicit ice generation?           : ', MP_doexpricit_icegen

    nstep_max = ceiling( ( TIME_DTSEC_ATMOS_PHY_MP * max_term_vel ) / minval( CDZ(:) ) )
    MP_ntmax_sedimentation = max( MP_ntmax_sedimentation, nstep_max )

    MP_NSTEP_SEDIMENTATION  = MP_ntmax_sedimentation
    MP_RNSTEP_SEDIMENTATION = 1.0_RP / real(MP_ntmax_sedimentation,kind=RP)
    MP_DTSEC_SEDIMENTATION  = TIME_DTSEC_ATMOS_PHY_MP * MP_RNSTEP_SEDIMENTATION

    if( IO_L ) write(IO_FID_LOG,*) '*** Timestep of sedimentation is divided into : ', MP_ntmax_sedimentation, 'step'
    if( IO_L ) write(IO_FID_LOG,*) '*** DT of sedimentation                       : ', MP_DTSEC_SEDIMENTATION, '[s]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_MP_TOMITA08,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_MP_TOMITA08. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_MP_TOMITA08)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** density of the snow    [kg/m3] : ', dens_s
    if( IO_L ) write(IO_FID_LOG,*) '*** density of the graupel [kg/m3] : ', dens_g
    if( IO_L ) write(IO_FID_LOG,*) '*** Nc for auto-conversion [num/m3]: ', autoconv_nc
    if( IO_L ) write(IO_FID_LOG,*) '*** Use k-k scheme?                : ', enable_kk2000
    if( IO_L ) write(IO_FID_LOG,*) '*** Use Roh scheme?                : ', enable_roh2014
    if( IO_L ) write(IO_FID_LOG,*)

    ATMOS_PHY_MP_tomita08_DENS(:)    = UNDEF
    ATMOS_PHY_MP_tomita08_DENS(I_HC) = dens_w
    ATMOS_PHY_MP_tomita08_DENS(I_HR) = dens_w
    ATMOS_PHY_MP_tomita08_DENS(I_HI) = dens_i
    ATMOS_PHY_MP_tomita08_DENS(I_HS) = dens_s
    ATMOS_PHY_MP_tomita08_DENS(I_HG) = dens_g

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

    if ( enable_roh2014 ) then ! overwrite parameters
       sw_roh2014 = 1.0_RP
       N0g        = 4.E8_RP
       As         = 0.069_RP
       Bs         = 2.0_RP
       Esi        = 0.25_RP
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

    do j = JS, JE
    do i = IS, IE
       Nc_def(i,j) = autoconv_nc
    enddo
    enddo

    if ( MP_doexpricit_icegen ) then
       only_liquid  = .true.
       sw_useicegen = 1.0_RP
    else
       only_liquid  = .false.
       sw_useicegen = 0.0_RP
    endif

    if ( enable_kk2000 ) then
       sw_kk2000 = 1.0_RP
    else
       sw_kk2000 = 0.0_RP
    endif

    return
  end subroutine ATMOS_PHY_MP_tomita08_setup

  !-----------------------------------------------------------------------------
  !> Cloud Microphysics
  subroutine ATMOS_PHY_MP_tomita08( &
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
    use scale_const, only: &
       DWATR => CONST_DWATR, &
       PI    => CONST_PI
    use scale_time, only: &
       dt => TIME_DTSEC_ATMOS_PHY_MP
    use scale_history, only: &
       HIST_in
    use scale_tracer, only: &
       QA, &
       TRACER_R, &
       TRACER_CV, &
       TRACER_MASS
    use scale_atmos_thermodyn, only: &
       THERMODYN_rhoe        => ATMOS_THERMODYN_rhoe,       &
       THERMODYN_rhot        => ATMOS_THERMODYN_rhot,       &
       THERMODYN_temp_pres_E => ATMOS_THERMODYN_temp_pres_E
    use scale_atmos_phy_mp_common, only: &
       MP_negative_fixer        => ATMOS_PHY_MP_negative_fixer,       &
       MP_precipitation         => ATMOS_PHY_MP_precipitation,        &
       MP_saturation_adjustment => ATMOS_PHY_MP_saturation_adjustment
    implicit none

    real(RP), intent(inout) :: DENS(KA,IA,JA)
    real(RP), intent(inout) :: MOMZ(KA,IA,JA)
    real(RP), intent(inout) :: MOMX(KA,IA,JA)
    real(RP), intent(inout) :: MOMY(KA,IA,JA)
    real(RP), intent(inout) :: RHOT(KA,IA,JA)
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA)
    real(RP), intent(in)    :: CCN(KA,IA,JA)
    real(RP), intent(out)   :: EVAPORATE(KA,IA,JA)   ! number of evaporated cloud [/m3/s]
    real(RP), intent(out)   :: SFLX_rain(IA,JA)
    real(RP), intent(out)   :: SFLX_snow(IA,JA)

    real(RP) :: RHOE_t(KA,IA,JA)
    real(RP) :: QTRC_t(KA,IA,JA,QA)
    real(RP) :: RHOE  (KA,IA,JA)
    real(RP) :: temp  (KA,IA,JA)
    real(RP) :: pres  (KA,IA,JA)

    real(RP) :: vterm     (KA,IA,JA,QA_MP-1) ! terminal velocity of each tracer [m/s]
    real(RP) :: FLX_rain  (KA,IA,JA)
    real(RP) :: FLX_snow  (KA,IA,JA)
    real(RP) :: wflux_rain(KA,IA,JA)
    real(RP) :: wflux_snow(KA,IA,JA)
    real(RP) :: qc_before_satadj(KA,IA,JA)
    integer  :: step
    integer  :: k, i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Cloud microphysics(tomita08)'

    if ( MP_donegative_fixer ) then
       call MP_negative_fixer( DENS(:,:,:),      & ! [INOUT]
                               RHOT(:,:,:),      & ! [INOUT]
                               QTRC(:,:,:,:),    & ! [INOUT]
                               I_QV,             & ! [IN]
                               MP_limit_negative ) ! [IN]
    endif

    call THERMODYN_rhoe( RHOE(:,:,:),   & ! [OUT]
                         RHOT(:,:,:),   & ! [IN]
                         QTRC(:,:,:,:), & ! [IN]
                         TRACER_CV(:),  & ! [IN]
                         TRACER_R(:),   & ! [IN]
                         TRACER_MASS(:) ) ! [IN]

    !##### MP Main #####
    call MP_tomita08( RHOE_t(:,:,:),   & ! [OUT]
                      QTRC_t(:,:,:,:), & ! [OUT]
                      RHOE  (:,:,:),   & ! [INOUT]
                      QTRC  (:,:,:,:), & ! [INOUT]
                      CCN   (:,:,:),   & ! [IN]
                      DENS  (:,:,:)    ) ! [IN]

    FLX_rain(:,:,:) = 0.0_RP
    FLX_snow(:,:,:) = 0.0_RP

    vterm(:,:,:,:) = 0.0_RP

    if ( MP_doprecipitation ) then

       do step = 1, MP_NSTEP_SEDIMENTATION

          call THERMODYN_temp_pres_E( temp(:,:,:),   & ! [OUT]
                                      pres(:,:,:),   & ! [OUT]
                                      DENS(:,:,:),   & ! [IN]
                                      RHOE(:,:,:),   & ! [IN]
                                      QTRC(:,:,:,:), & ! [IN]
                                      TRACER_CV(:),  & ! [IN]
                                      TRACER_R(:),   & ! [IN]
                                      TRACER_MASS(:) ) ! [IN]

          call MP_tomita08_vterm( vterm(:,:,:,:), & ! [OUT]
                                  DENS (:,:,:),   & ! [IN]
                                  temp (:,:,:),   & ! [IN]
                                  QTRC (:,:,:,:)  ) ! [IN]

          call MP_precipitation( wflux_rain(:,:,:),     & ! [OUT]
                                 wflux_snow(:,:,:),     & ! [OUT]
                                 DENS    (:,:,:),       & ! [INOUT]
                                 MOMZ    (:,:,:),       & ! [INOUT]
                                 MOMX    (:,:,:),       & ! [INOUT]
                                 MOMY    (:,:,:),       & ! [INOUT]
                                 RHOE    (:,:,:),       & ! [INOUT]
                                 QTRC    (:,:,:,:),     & ! [INOUT]
                                 QA_MP,                 & ! [IN]
                                 QS_MP,                 & ! [IN]
                                 vterm   (:,:,:,:),     & ! [IN]
                                 temp    (:,:,:),       & ! [IN]
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

    else
       vterm(:,:,:,:) = 0.0_RP
    endif

    call HIST_in( vterm(:,:,:,I_mp_QR), 'Vterm_QR', 'terminal velocity of QR', 'm/s' )
    call HIST_in( vterm(:,:,:,I_mp_QI), 'Vterm_QI', 'terminal velocity of QI', 'm/s' )
    call HIST_in( vterm(:,:,:,I_mp_QS), 'Vterm_QS', 'terminal velocity of QS', 'm/s' )
    call HIST_in( vterm(:,:,:,I_mp_QG), 'Vterm_QG', 'terminal velocity of QG', 'm/s' )

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       qc_before_satadj(k,i,j) = QTRC(k,i,j,I_QC)
    enddo
    enddo
    enddo

    call MP_saturation_adjustment( RHOE_t(:,:,:),   & ! [INOUT]
                                   QTRC_t(:,:,:,:), & ! [INOUT]
                                   RHOE  (:,:,:),   & ! [INOUT]
                                   QTRC  (:,:,:,:), & ! [INOUT]
                                   DENS  (:,:,:),   & ! [IN]
                                   I_QV,            & ! [IN]
                                   I_QC,            & ! [IN]
                                   I_QI,            & ! [IN]
                                   only_liquid      ) ! [IN]

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       EVAPORATE(k,i,j) = max( qc_before_satadj(k,i,j) - QTRC(k,i,j,I_QC), 0.0_RP ) / dt ! if negative, condensation
       EVAPORATE(k,i,j) = EVAPORATE(k,i,j) * DENS(k,i,j) &
                        / (4.0_RP/3.0_RP*PI*DWATR*re_qc**3)  ! mass -> number (assuming constant particle radius as re_qc)
    enddo
    enddo
    enddo

    !##### END MP Main #####

    call THERMODYN_rhot( RHOT(:,:,:),   & ! [OUT]
                         RHOE(:,:,:),   & ! [IN]
                         QTRC(:,:,:,:), & ! [IN]
                         TRACER_CV(:),  & ! [IN]
                         TRACER_R(:),   & ! [IN]
                         TRACER_MASS(:) ) ! [IN]

    if ( MP_donegative_fixer ) then
       call MP_negative_fixer( DENS(:,:,:),      & ! [INOUT]
                               RHOT(:,:,:),      & ! [INOUT]
                               QTRC(:,:,:,:),    & ! [INOUT]
                               I_QV,             & ! [IN]
                               MP_limit_negative ) ! [IN]
    endif

    SFLX_rain(:,:) = FLX_rain(KS-1,:,:)
    SFLX_snow(:,:) = FLX_snow(KS-1,:,:)

    return
  end subroutine ATMOS_PHY_MP_tomita08

  !-----------------------------------------------------------------------------
  !> Lin-type cold rain microphysics
  subroutine MP_tomita08( &
       RHOE_t,    &
       QTRC_t,    &
       RHOE0,     &
       QTRC0,     &
       CCN,       &
       DENS0      )
    use scale_const, only: &
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
    use scale_time, only: &
       dt => TIME_DTSEC_ATMOS_PHY_MP
    use scale_tracer, only: &
       QA, &
       TRACER_R, &
       TRACER_CV, &
       TRACER_MASS
    use scale_history, only: &
       HIST_in
    use scale_atmos_thermodyn, only: &
       THERMODYN_temp_pres_E => ATMOS_THERMODYN_temp_pres_E
    use scale_atmos_hydrometeor, only: &
       LHV, &
       LHF
    use scale_atmos_saturation, only: &
       SATURATION_dens2qsat_liq => ATMOS_SATURATION_dens2qsat_liq, &
       SATURATION_dens2qsat_ice => ATMOS_SATURATION_dens2qsat_ice
    implicit none

    real(RP), intent(out)   :: RHOE_t(KA,IA,JA)    ! tendency rhoe             [J/m3/s]
    real(RP), intent(out)   :: QTRC_t(KA,IA,JA,QA) ! tendency tracer           [kg/kg/s]
    real(RP), intent(inout) :: RHOE0 (KA,IA,JA)    ! density * internal energy [J/m3]
    real(RP), intent(inout) :: QTRC0 (KA,IA,JA,QA) ! mass concentration        [kg/kg]
    real(RP), intent(in)    :: CCN   (KA,IA,JA)    ! CCN number concentration  [#/m3]
    real(RP), intent(in)    :: DENS0 (KA,IA,JA)    ! density                   [kg/m3]

    ! working
    real(RP) :: rdt

    real(RP) :: TEMP0(KA,IA,JA) ! temperature [K]
    real(RP) :: PRES0(KA,IA,JA) ! pressure    [Pa]
    real(RP) :: QSATL(KA,IA,JA) ! saturated water vapor for liquid water [kg/kg]
    real(RP) :: QSATI(KA,IA,JA) ! saturated water vapor for ice water    [kg/kg]

    real(RP) :: dens
    real(RP) :: rhoe
    real(RP) :: temp
    real(RP) :: pres
    real(RP) :: q(QS_MP:QE_MP)
    real(RP) :: Sliq  ! saturated ratio S for liquid water [0-1]
    real(RP) :: Sice  ! saturated ratio S for ice water    [0-1]

    real(RP) :: Rdens
    real(RP) :: rho_fact ! density factor
    real(RP) :: temc     ! T - T0 [K]

    real(RP) :: RLMDr, RLMDr_2, RLMDr_3
    real(RP) :: RLMDs, RLMDs_2, RLMDs_3
    real(RP) :: RLMDg, RLMDg_2, RLMDg_3
    real(RP) :: RLMDr_1br, RLMDr_2br, RLMDr_3br
    real(RP) :: RLMDs_1bs, RLMDs_2bs, RLMDs_3bs
    real(RP) :: RLMDr_dr, RLMDr_3dr, RLMDr_5dr
    real(RP) :: RLMDs_ds, RLMDs_3ds, RLMDs_5ds
    real(RP) :: RLMDg_dg, RLMDg_3dg, RLMDg_5dg
    real(RP) :: RLMDr_7
    real(RP) :: RLMDr_6dr

    !---< Roh and Satoh (2014) >---
    real(RP) :: tems, rhoqs, Xs2
    real(RP) :: MOMs_0, MOMs_1, MOMs_2
    real(RP) :: MOMs_0bs, MOMs_1bs, MOMs_2bs
    real(RP) :: MOMs_2ds, MOMs_5ds_h, RMOMs_Vt
    real(RP) :: coef_at(4), coef_bt(4)
    real(RP) :: loga_, b_, nm

    real(RP) :: Vti, Vtr, Vts, Vtg  !< terminal velocity
    real(RP) :: Esi_mod, Egs_mod    !< modified accretion efficiency
    real(RP) :: rhoqc               !< rho * qc
    real(RP) :: Nc(KA,IA,JA)        !< Number concentration of cloud water [1/cc]
    real(RP) :: Pracw_orig,  Pracw_kk !< accretion       term by orig  & k-k scheme
    real(RP) :: Praut_berry, Praut_kk !< auto-conversion term by berry & k-k scheme
    real(RP) :: Dc                  !< relative variance
    real(RP) :: betai, betas        !< sticky parameter for auto-conversion
    real(RP) :: Da                  !< thermal diffusion coefficient of air
    real(RP) :: Kd                  !< diffusion coefficient of water vapor in air
    real(RP) :: Nu                  !< kinematic viscosity of air
    real(RP) :: Glv, Giv, Gil       !< thermodynamic function
    real(RP) :: ventr, vents, ventg !< ventilation factor
    real(RP) :: Bergeron_sw         !< if 0C<T<30C, sw=1
    real(RP) :: a1 (KA,IA,JA)       !<
    real(RP) :: a2 (KA,IA,JA)       !<
    real(RP) :: ma2(KA,IA,JA)       !< 1-a2
    real(RP) :: dt1                 !< time during which the an ice particle of 40um grows to 50um
    real(RP) :: Ni50                !< number concentration of ice particle of 50um

    ! for explicit ice generation
    real(RP), parameter :: Di_max = 500.E-6_RP
    real(RP), parameter :: Di_a   = 11.9_RP
    real(RP) :: sw, rhoqi, XNi, XMi, Di, Ni0, Qi0, dq, qtmp
    real(RP) :: Pidep, Pigen

    real(RP) :: ice

    real(RP) :: net, fac, fac_sw, fac_warm, fac_ice
    real(RP) :: w(w_nmax)

    real(RP) :: tend(QS_MP:QE_MP)

    real(RP) :: zerosw
    real(RP) :: tmp

    integer  :: k, i, j, iq
    !---------------------------------------------------------------------------

    call PROF_rapstart('MP_tomita08', 3)

!     RHOE_t(:,:,:)   = 0.0_RP
!     QTRC_t(:,:,:,:) = 0.0_RP

    if ( MP_couple_aerosol ) then
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          Nc(k,i,j) = max( CCN(k,i,j)*1.E-6_RP, Nc_def(i,j) ) ! [#/m3]->[#/cc]
       enddo
       enddo
       enddo
    else
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          Nc(k,i,j) = Nc_def(i,j)
       enddo
       enddo
       enddo
    endif

    rdt = 1.0_RP / dt

    call THERMODYN_temp_pres_E( TEMP0(:,:,:),   & ! [OUT]
                                PRES0(:,:,:),   & ! [OUT]
                                DENS0(:,:,:),   & ! [IN]
                                RHOE0(:,:,:),   & ! [IN]
                                QTRC0(:,:,:,:), & ! [IN]
                                TRACER_CV(:),   & ! [IN]
                                TRACER_R(:),    & ! [IN]
                                TRACER_MASS(:)  ) ! [IN]

    call SATURATION_dens2qsat_liq( QSATL(:,:,:), & ! [OUT]
                                   TEMP0(:,:,:), & ! [IN]
                                   DENS0(:,:,:)  ) ! [IN]

    call SATURATION_dens2qsat_ice( QSATI(:,:,:), & ! [OUT]
                                   TEMP0(:,:,:), & ! [IN]
                                   DENS0(:,:,:)  ) ! [IN]

    call MP_tomita08_BergeronParam( TEMP0(:,:,:), & ! [IN]
                                    a1   (:,:,:), & ! [OUT]
                                    a2   (:,:,:), & ! [OUT]
                                    ma2  (:,:,:)  ) ! [OUT]

!$omp parallel do &
!$omp private(tend, coef_bt, coef_at, q, w) &
!$omp private(dens, rhoe, temp, pres) &
!$omp private(Sliq, Sice, Rdens, rho_fact, temc) &
!$omp private(Bergeron_sw, zerosw) &
!$omp private(RLMDr, RLMDr_dr, RLMDr_2, RLMDr_3, RLMDr_7, RLMDr_1br, RLMDr_2br, RLMDr_3br, RLMDr_3dr, RLMDr_5dr, RLMDr_6dr) &
!$omp private(RLMDs, RLMDs_ds, RLMDs_2, RLMDs_3,          RLMDs_1bs, RLMDs_2bs, RLMDs_3bs, RLMDs_3ds, RLMDs_5ds) &
!$omp private(MOMs_0, MOMs_1, MOMs_2, MOMs_0bs, MOMs_1bs, MOMs_2bs, MOMs_2ds, MOMs_5ds_h, RMOMs_Vt) &
!$omp private(rhoqc, rhoqs, Xs2, tems, loga_, b_, nm) &
!$omp private(RLMDg, RLMDg_dg, RLMDg_2, RLMDg_3, RLMDg_3dg, RLMDg_5dg) &
!$omp private(Vti, Vtr, Vts, Vtg, Esi_mod, Egs_mod, Dc, Praut_berry, Praut_kk, betai, betas, Da, Kd, NU, Glv, Giv, Gil, ventr, vents, ventg) &
!$omp private(tmp, dt1, Ni50, ice, net, fac_sw, fac, fac_warm, fac_ice) &
!$omp collapse(3)
!OCL TEMP_PRIVATE(tend, coef_bt, coef_at, q, w)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE

       ! store to work
       dens = DENS0(k,i,j)
       rhoe = RHOE0(k,i,j)
       temp = TEMP0(k,i,j)
       pres = PRES0(k,i,j)
       do iq = QS_MP, QE_MP
          q(iq) = QTRC0(k,i,j,iq)
       enddo

       ! saturation ratio S
       Sliq = q(I_QV) / max( QSATL(k,i,j), EPS )
       Sice = q(I_QV) / max( QSATI(k,i,j), EPS )

       Rdens    = 1.0_RP / dens
       rho_fact = sqrt( dens00 * Rdens )
       temc     = temp - TEM00

       w(I_delta1) = ( 0.5_RP + sign(0.5_RP, q(I_QR) - 1.E-4_RP ) )

       w(I_delta2) = ( 0.5_RP + sign(0.5_RP, 1.E-4_RP - q(I_QR) ) ) &
                   * ( 0.5_RP + sign(0.5_RP, 1.E-4_RP - q(I_QS) ) )

       w(I_delta3) = 0.5_RP + sign(0.5_RP, Sice - 1.0_RP )

       w(I_dqv_dt) = q(I_QV) * rdt
       w(I_dqc_dt) = q(I_QC) * rdt
       w(I_dqr_dt) = q(I_QR) * rdt
       w(I_dqi_dt) = q(I_QI) * rdt
       w(I_dqs_dt) = q(I_QS) * rdt
       w(I_dqg_dt) = q(I_QG) * rdt

       Bergeron_sw = ( 0.5_RP + sign(0.5_RP, temc + 30.0_RP ) ) &
                   * ( 0.5_RP + sign(0.5_RP, 0.0_RP - temc  ) ) &
                   * ( 1.0_RP - sw_useicegen                  )

       ! slope parameter lambda (Rain)
       zerosw = 0.5_RP - sign(0.5_RP, q(I_QR) - 1.E-12_RP )
       RLMDr  = sqrt(sqrt( dens * q(I_QR) / ( Ar * N0r * GAM_1br ) + zerosw )) * ( 1.0_RP - zerosw )

       RLMDr_dr  = sqrt( RLMDr )       ! **Dr
       RLMDr_2   = RLMDr**2
       RLMDr_3   = RLMDr**3
       RLMDr_7   = RLMDr**7
       RLMDr_1br = RLMDr**4 ! (1+Br)
       RLMDr_2br = RLMDr**5 ! (2+Br)
       RLMDr_3br = RLMDr**6 ! (3+Br)
       RLMDr_3dr = RLMDr**3 * RLMDr_dr
       RLMDr_5dr = RLMDr**5 * RLMDr_dr
       RLMDr_6dr = RLMDr**6 * RLMDr_dr

       ! slope parameter lambda (Snow)
       zerosw = 0.5_RP - sign(0.5_RP, q(I_QS) - 1.E-12_RP )
       RLMDs  = sqrt(sqrt( dens * q(I_QS) / ( As * N0s * GAM_1bs ) + zerosw )) * ( 1.0_RP - zerosw )

       RLMDs_ds  = sqrt( sqrt(RLMDs) ) ! **Ds
       RLMDs_2   = RLMDs**2
       RLMDs_3   = RLMDs**3
       RLMDs_1bs = RLMDs**4 ! (1+Bs)
       RLMDs_2bs = RLMDs**5 ! (2+Bs)
       RLMDs_3bs = RLMDs**6 ! (3+Bs)
       RLMDs_3ds = RLMDs**3 * RLMDs_ds
       RLMDs_5ds = RLMDs**5 * RLMDs_ds

       MOMs_0     = N0s * GAM       * RLMDs           ! Ns * 0th moment
       MOMs_1     = N0s * GAM_2     * RLMDs_2         ! Ns * 1st moment
       MOMs_2     = N0s * GAM_3     * RLMDs_3         ! Ns * 2nd moment
       MOMs_0bs   = N0s * GAM_1bs   * RLMDs_1bs       ! Ns * 0+bs
       MOMs_1bs   = N0s * GAM_2bs   * RLMDs_2bs       ! Ns * 1+bs
       MOMs_2bs   = N0s * GAM_3bs   * RLMDs_3bs       ! Ns * 2+bs
       MOMs_2ds   = N0s * GAM_3ds   * RLMDs_3ds       ! Ns * 2+ds
       MOMs_5ds_h = N0s * GAM_5ds_h * sqrt(RLMDs_5ds) ! Ns * (5+ds)/2
       RMOMs_Vt   = GAM_1bsds / GAM_1bs * RLMDs_ds

       !---< modification by Roh and Satoh (2014) >---
       ! bimodal size distribution of snow
       rhoqs  = dens * q(I_QS)
       zerosw = 0.5_RP - sign(0.5_RP, rhoqs - 1.E-12_RP )
       Xs2    = rhoqs / As * ( 1.0_RP - zerosw )

       tems       = min( -0.1_RP, temc )
       coef_at(1) = coef_a( 1) + tems * ( coef_a( 2) + tems * ( coef_a( 5) + tems * coef_a( 9) ) )
       coef_at(2) = coef_a( 3) + tems * ( coef_a( 4) + tems *   coef_a( 7) )
       coef_at(3) = coef_a( 6) + tems *   coef_a( 8)
       coef_at(4) = coef_a(10)
       coef_bt(1) = coef_b( 1) + tems * ( coef_b( 2) + tems * ( coef_b( 5) + tems * coef_b( 9) ) )
       coef_bt(2) = coef_b( 3) + tems * ( coef_b( 4) + tems *   coef_b( 7) )
       coef_bt(3) = coef_b( 6) + tems *   coef_b( 8)
       coef_bt(4) = coef_b(10)
       ! 0th moment
       loga_  = coef_at(1)
       b_     = coef_bt(1)
       MOMs_0 = (        sw_roh2014 ) * 10.0_RP**loga_ * Xs2**b_ &
              + ( 1.0_RP-sw_roh2014 ) * MOMs_0
       ! 1st moment
       nm = 1.0_RP
       loga_  = coef_at(1) + nm * ( coef_at(2) + nm * ( coef_at(3) + nm * coef_at(4) ) )
          b_  = coef_bt(1) + nm * ( coef_bt(2) + nm * ( coef_bt(3) + nm * coef_bt(4) ) )
       MOMs_1 = (        sw_roh2014 ) * 10.0_RP**loga_ * Xs2**b_ &
              + ( 1.0_RP-sw_roh2014 ) * MOMs_1
       ! 2nd moment
       MOMs_2 = (        sw_roh2014 ) * Xs2 &
              + ( 1.0_RP-sw_roh2014 ) * MOMs_2
       ! 0 + Bs(=2) moment
       nm = 2.0_RP
       loga_  = coef_at(1) + nm * ( coef_at(2) + nm * ( coef_at(3) + nm * coef_at(4) ) )
          b_  = coef_bt(1) + nm * ( coef_bt(2) + nm * ( coef_bt(3) + nm * coef_bt(4) ) )
       MOMs_0bs = (        sw_roh2014 ) * 10.0_RP**loga_ * Xs2**b_ &
                + ( 1.0_RP-sw_roh2014 ) * MOMs_0bs
       ! 1 + Bs(=2) moment
       nm = 3.0_RP
       loga_  = coef_at(1) + nm * ( coef_at(2) + nm * ( coef_at(3) + nm * coef_at(4) ) )
          b_  = coef_bt(1) + nm * ( coef_bt(2) + nm * ( coef_bt(3) + nm * coef_bt(4) ) )
       MOMs_1bs = (        sw_roh2014 ) * 10.0_RP**loga_ * Xs2**b_ &
                + ( 1.0_RP-sw_roh2014 ) * MOMs_1bs
       ! 2 + Bs(=2) moment
       nm = 4.0_RP
       loga_  = coef_at(1) + nm * ( coef_at(2) + nm * ( coef_at(3) + nm * coef_at(4) ) )
          b_  = coef_bt(1) + nm * ( coef_bt(2) + nm * ( coef_bt(3) + nm * coef_bt(4) ) )
       MOMs_2bs = (        sw_roh2014 ) * 10.0_RP**loga_ * Xs2**b_ &
                + ( 1.0_RP-sw_roh2014 ) * MOMs_2bs
       ! 2 + Ds(=0.25) moment
       nm = 2.25_RP
       loga_  = coef_at(1) + nm * ( coef_at(2) + nm * ( coef_at(3) + nm * coef_at(4) ) )
          b_  = coef_bt(1) + nm * ( coef_bt(2) + nm * ( coef_bt(3) + nm * coef_bt(4) ) )
       MOMs_2ds = (        sw_roh2014 ) * 10.0_RP**loga_ * Xs2**b_ &
                + ( 1.0_RP-sw_roh2014 ) * MOMs_2ds
       ! ( 3 + Ds(=0.25) ) / 2  moment
       nm = 1.625_RP
       loga_  = coef_at(1) + nm * ( coef_at(2) + nm * ( coef_at(3) + nm * coef_at(4) ) )
          b_  = coef_bt(1) + nm * ( coef_bt(2) + nm * ( coef_bt(3) + nm * coef_bt(4) ) )
       MOMs_5ds_h = (        sw_roh2014 ) * 10.0_RP**loga_ * Xs2**b_ &
                  + ( 1.0_RP-sw_roh2014 ) * MOMs_5ds_h
       ! Bs(=2) + Ds(=0.25) moment
       nm = 2.25_RP
       loga_  = coef_at(1) + nm * ( coef_at(2) + nm * ( coef_at(3) + nm * coef_at(4) ) )
          b_  = coef_bt(1) + nm * ( coef_bt(2) + nm * ( coef_bt(3) + nm * coef_bt(4) ) )
       RMOMs_Vt = (        sw_roh2014 ) * 10.0_RP**loga_ * Xs2**b_ / ( MOMs_0bs + zerosw ) &
                + ( 1.0_RP-sw_roh2014 ) * RMOMs_Vt

       ! slope parameter lambda (Graupel)
       zerosw = 0.5_RP - sign(0.5_RP, q(I_QG) - 1.E-12_RP )
       RLMDg  = sqrt(sqrt( dens * q(I_QG) / ( Ag * N0g * GAM_1bg ) + zerosw )) * ( 1.0_RP - zerosw )

       RLMDg_dg  = sqrt( RLMDg )       ! **Dg
       RLMDg_2   = RLMDg**2
       RLMDg_3   = RLMDg**3
       RLMDg_3dg = RLMDg**3 * RLMDg_dg
       RLMDg_5dg = RLMDg**5 * RLMDg_dg

       w(I_RLMDr) = RLMDr
       w(I_RLMDs) = RLMDs
       w(I_RLMDg) = RLMDg

       !---< terminal velocity >
       zerosw = 0.5_RP + sign(0.5_RP, q(I_QI) - 1.E-8_RP )
       Vti = -3.29_RP * ( dens * q(I_QI) * zerosw )**0.16_RP
       Vtr = -Cr * rho_fact * GAM_1brdr / GAM_1br * RLMDr_dr
       Vts = -Cs * rho_fact * RMOMs_Vt
       Vtg = -Cg * rho_fact * GAM_1bgdg / GAM_1bg * RLMDg_dg

       !---< Accretion >---
       Esi_mod = min( Esi, Esi * exp( gamma_sacr * temc ) )
       Egs_mod = min( Egs, Egs * exp( gamma_gacs * temc ) )

       ! [Pracw] accretion rate of cloud water by rain
       Pracw_orig = q(I_QC) * 0.25_RP * PI * Erw * N0r * Cr * GAM_3dr * RLMDr_3dr * rho_fact

       Pracw_kk   = 67.0_RP * ( q(I_QC) * q(I_QR) )**1.15_RP ! eq.(33) in KK(2000)

       ! switch orig / k-k scheme
       w(I_Pracw) = ( 1.0_RP - sw_kk2000 ) * Pracw_orig &
                  + (          sw_kk2000 ) * Pracw_kk

       ! [Psacw] accretion rate of cloud water by snow
       w(I_Psacw) = q(I_QC) * 0.25_RP * PI * Esw       * Cs * MOMs_2ds            * rho_fact

       ! [Pgacw] accretion rate of cloud water by graupel
       w(I_Pgacw) = q(I_QC) * 0.25_RP * PI * Egw * N0g * Cg * GAM_3dg * RLMDg_3dg * rho_fact

       ! [Praci] accretion rate of cloud ice by rain
       w(I_Praci) = q(I_QI) * 0.25_RP * PI * Eri * N0r * Cr * GAM_3dr * RLMDr_3dr * rho_fact

       ! [Psaci] accretion rate of cloud ice by snow
       w(I_Psaci) = q(I_QI) * 0.25_RP * PI * Esi_mod   * Cs * MOMs_2ds            * rho_fact

       ! [Pgaci] accretion rate of cloud ice by grupel
       w(I_Pgaci) = q(I_QI) * 0.25_RP * PI * Egi * N0g * Cg * GAM_3dg * RLMDg_3dg * rho_fact * ( 1.0_RP-sw_roh2014 )

       ! [Piacr] accretion rate of rain by cloud ice
       w(I_Piacr) = q(I_QI) * Ar / 4.19E-13_RP * 0.25_RP * PI * Eri * N0r * Cr * GAM_6dr * RLMDr_6dr * rho_fact

       ! [Psacr] accretion rate of rain by snow
       w(I_Psacr) = Ar * 0.25_RP * PI * Rdens * Esr * N0r       * abs(Vtr-Vts) &
                  * (          GAM_1br * RLMDr_1br * MOMs_2          &
                    + 2.0_RP * GAM_2br * RLMDr_2br * MOMs_1          &
                    +          GAM_3br * RLMDr_3br * MOMs_0          )

       ! [Pgacr] accretion rate of rain by graupel
       w(I_Pgacr) = Ar * 0.25_RP * PI * Rdens * Egr * N0g * N0r * abs(Vtg-Vtr) &
                  * (          GAM_1br * RLMDr_1br * GAM_3 * RLMDg_3 &
                    + 2.0_RP * GAM_2br * RLMDr_2br * GAM_2 * RLMDg_2 &
                    +          GAM_3br * RLMDr_3br * GAM   * RLMDg   )

       ! [Pracs] accretion rate of snow by rain
       w(I_Pracs) = As * 0.25_RP * PI * Rdens * Esr       *  N0r * abs(Vtr-Vts) &
                  * (          MOMs_0bs            * GAM_3 * RLMDr_3 &
                    + 2.0_RP * MOMs_1bs            * GAM_2 * RLMDr_2 &
                    +          MOMs_2bs            * GAM   * RLMDr   )

       ! [Pgacs] accretion rate of snow by graupel
       w(I_Pgacs) = As * 0.25_RP * PI * Rdens * Egs_mod   * N0g * abs(Vtg-Vts) * ( 1.0_RP-sw_roh2014 ) &
                  * (          MOMs_0bs            * GAM_3 * RLMDg_3 &
                    + 2.0_RP * MOMs_1bs            * GAM_2 * RLMDg_2 &
                    +          MOMs_2bs            * GAM   * RLMDg   )

       !---< Auto-conversion >---
       ! [Praut] auto-conversion rate from cloud water to rain
       rhoqc = dens * q(I_QC) * 1000.0_RP ! [g/m3]
       Dc    = 0.146_RP - 5.964E-2_RP * log( Nc(k,i,j) / 2000.0_RP )
       Praut_berry = Rdens * 1.67E-5_RP * rhoqc * rhoqc / ( 5.0_RP + 3.66E-2_RP * Nc(k,i,j) / ( Dc * rhoqc + EPS ) )

       Praut_kk    = 1350.0_RP * q(I_QC)**2.47_RP * Nc(k,i,j)**(-1.79_RP) ! eq.(29) in KK(2000)

       ! switch berry / k-k scheme
       w(I_Praut) = ( 1.0_RP - sw_kk2000 ) * Praut_berry &
                  + (          sw_kk2000 ) * Praut_kk

       ! [Psaut] auto-conversion rate from cloud ice to snow
       betai = min( beta_saut, beta_saut * exp( gamma_saut * temc ) )
       w(I_Psaut) = max( betai*(q(I_QI)-qicrt_saut), 0.0_RP )

       ! [Pgaut] auto-conversion rate from snow to graupel
       betas = min( beta_gaut, beta_gaut * exp( gamma_gaut * temc ) )
       w(I_Pgaut) = max( betas*(q(I_QS)-qscrt_gaut), 0.0_RP )

       !---< Evaporation, Sublimation >---
       Da = ( Da0 + dDa_dT * temc )
       Kd = ( Dw0 + dDw_dT * temc ) * PRE00 / pres
       NU = ( mu0 + dmu_dT * temc ) * Rdens

       Glv  = 1.0_RP / ( LHV0/(Da*temp) * ( LHV0/(Rvap*temp) - 1.0_RP ) + 1.0_RP/(Kd*dens*QSATL(k,i,j)) )
       Giv  = 1.0_RP / ( LHS0/(Da*temp) * ( LHS0/(Rvap*temp) - 1.0_RP ) + 1.0_RP/(Kd*dens*QSATI(k,i,j)) )
       Gil  = 1.0_RP / LHF0 * (Da*temc)

       ! [Prevp] evaporation rate of rain
       ventr = f1r * GAM_2 * RLMDr_2 + f2r * sqrt( Cr * rho_fact / NU * RLMDr_5dr ) * GAM_5dr_h

       w(I_Prevp) = 2.0_RP * PI * Rdens * N0r * ( 1.0_RP-min(Sliq,1.0_RP) ) * Glv * ventr

       ! [Psdep,Pssub] deposition/sublimation rate for snow
       vents = f1s * MOMs_1          + f2s * sqrt( Cs * rho_fact / NU             ) * MOMs_5ds_h

       tmp = 2.0_RP * PI * Rdens *       ( Sice-1.0_RP ) * Giv * vents

       w(I_Psdep) = ( w(I_delta3)        ) * tmp ! Sice < 1
       w(I_Pssub) = ( w(I_delta3)-1.0_RP ) * tmp ! Sice > 1

       ! [Psmlt] melting rate of snow
       w(I_Psmlt) = 2.0_RP * PI * Rdens *       Gil * vents &
                  + CL * temc / LHF0 * ( w(I_Psacw) + w(I_Psacr) )

       ! [Pgdep/pgsub] deposition/sublimation rate for graupel
       ventg = f1g * GAM_2 * RLMDg_2 + f2g * sqrt( Cg * rho_fact / NU * RLMDg_5dg ) * GAM_5dg_h

       tmp = 2.0_RP * PI * Rdens * N0g * ( Sice-1.0_RP ) * Giv * ventg

       w(I_Pgdep) = ( w(I_delta3)        ) * tmp ! Sice < 1
       w(I_Pgsub) = ( w(I_delta3)-1.0_RP ) * tmp ! Sice > 1

       ! [Pgmlt] melting rate of graupel
       w(I_Pgmlt) = 2.0_RP * PI * Rdens * N0g * Gil * ventg &
                  + CL * temc / LHF0 * ( w(I_Pgacw) + w(I_Pgacr) )

       ! [Pgfrz] freezing rate of graupel
       w(I_Pgfrz) = 2.0_RP * PI * Rdens * N0r * 60.0_RP * B_frz * Ar * ( exp(-A_frz*temc) - 1.0_RP ) * RLMDr_7

       ! [Psfw,Psfi] ( Bergeron process ) growth rate of snow by Bergeron process from cloud water/ice
       dt1  = ( mi50**ma2(k,i,j) - mi40**ma2(k,i,j) ) / ( a1(k,i,j) * ma2(k,i,j) )
       Ni50 = q(I_QI) * dt / ( mi50 * dt1 )

       w(I_Psfw) = Bergeron_sw * Ni50 * ( a1(k,i,j) * mi50**a2(k,i,j)                   &
                                        + PI * Eiw * dens * q(I_QC) * Ri50*Ri50 * vti50 )
       w(I_Psfi) = Bergeron_sw * q(I_QI) / dt1

       w(I_Psdep) = min( w(I_Psdep), w(I_dqv_dt) )
       w(I_Pgdep) = min( w(I_Pgdep), w(I_dqv_dt) )

       w(I_Praut) = min( w(I_Praut), w(I_dqc_dt) )
       w(I_Pracw) = min( w(I_Pracw), w(I_dqc_dt) )
       w(I_Psacw) = min( w(I_Psacw), w(I_dqc_dt) )
       w(I_Pgacw) = min( w(I_Pgacw), w(I_dqc_dt) )
       w(I_Psfw ) = min( w(I_Psfw ), w(I_dqc_dt) )

       w(I_Prevp) = min( w(I_Prevp), w(I_dqr_dt) )
       w(I_Piacr) = min( w(I_Piacr), w(I_dqr_dt) )
       w(I_Psacr) = min( w(I_Psacr), w(I_dqr_dt) )
       w(I_Pgacr) = min( w(I_Pgacr), w(I_dqr_dt) )
       w(I_Pgfrz) = min( w(I_Pgfrz), w(I_dqr_dt) )

       w(I_Psaut) = min( w(I_Psaut), w(I_dqi_dt) )
       w(I_Praci) = min( w(I_Praci), w(I_dqi_dt) )
       w(I_Psaci) = min( w(I_Psaci), w(I_dqi_dt) )
       w(I_Pgaci) = min( w(I_Pgaci), w(I_dqi_dt) )
       w(I_Psfi ) = min( w(I_Psfi ), w(I_dqi_dt) )

       w(I_Pgaut) = min( w(I_Pgaut), w(I_dqs_dt) )
       w(I_Pracs) = min( w(I_Pracs), w(I_dqs_dt) )
       w(I_Pgacs) = min( w(I_Pgacs), w(I_dqs_dt) )
       w(I_Psmlt) = max( w(I_Psmlt), 0.0_RP      )
       w(I_Psmlt) = min( w(I_Psmlt), w(I_dqs_dt) )
       w(I_Pssub) = min( w(I_Pssub), w(I_dqs_dt) )

       w(I_Pgmlt) = max( w(I_Pgmlt), 0.0_RP        )
       w(I_Pgmlt) = min( w(I_Pgmlt), w(I_dqg_dt) )
       w(I_Pgsub) = min( w(I_Pgsub), w(I_dqg_dt) )

       w(I_Piacr_s) = ( 1.0_RP - w(I_delta1) ) * w(I_Piacr)
       w(I_Piacr_g) = (          w(I_delta1) ) * w(I_Piacr)
       w(I_Praci_s) = ( 1.0_RP - w(I_delta1) ) * w(I_Praci)
       w(I_Praci_g) = (          w(I_delta1) ) * w(I_Praci)
       w(I_Psacr_s) = (          w(I_delta2) ) * w(I_Psacr)
       w(I_Psacr_g) = ( 1.0_RP - w(I_delta2) ) * w(I_Psacr)
       w(I_Pracs  ) = ( 1.0_RP - w(I_delta2) ) * w(I_Pracs)

       ice = 0.5_RP - sign( 0.5_RP, TEMP0(k,i,j)-TEM00 ) ! 0: warm, 1: ice

       ! [QC]
       net = &
           - w(I_Praut) & ! [loss] c->r
           - w(I_Pracw) & ! [loss] c->r
           - w(I_Psacw) & ! [loss] c->s
           - w(I_Pgacw) & ! [loss] c->g
           - w(I_Psfw ) * ice ! [loss] c->s

       fac_sw = 0.5_RP + sign( 0.5_RP, net+EPS ) ! if production > loss , fac_sw=1
       fac    = (          fac_sw ) &
              + ( 1.0_RP - fac_sw ) * min( -w(I_dqc_dt)/(net-fac_sw), 1.0_RP ) ! loss limiter

       w(I_Praut) = w(I_Praut) * fac
       w(I_Pracw) = w(I_Pracw) * fac
       w(I_Psacw) = w(I_Psacw) * fac
       w(I_Pgacw) = w(I_Pgacw) * fac
       w(I_Psfw ) = w(I_Psfw ) * fac

       ! [QI]
       net = &
           - w(I_Psaut  ) & ! [loss] i->s
           - w(I_Praci_s) & ! [loss] i->s
           - w(I_Psaci  ) & ! [loss] i->s
           - w(I_Psfi   ) & ! [loss] i->s
           - w(I_Praci_g) & ! [loss] i->g
           - w(I_Pgaci  )   ! [loss] i->g

       fac_sw = 0.5_RP + sign( 0.5_RP, net+EPS ) ! if production > loss , fac_sw=1
       fac    = (          fac_sw ) &
              + ( 1.0_RP - fac_sw ) * min( -w(I_dqi_dt)/(net-fac_sw), 1.0_RP ) ! loss limiter

       w(I_Praci_s) = w(I_Praci_s) * fac
       w(I_Psaut  ) = w(I_Psaut  ) * fac
       w(I_Psaci  ) = w(I_Psaci  ) * fac
       w(I_Psfi   ) = w(I_Psfi   ) * fac
       w(I_Praci_g) = w(I_Praci_g) * fac
       w(I_Pgaci  ) = w(I_Pgaci  ) * fac

       ! [QR]
       net = w(I_Praut  ) & ! [prod] c->r
           + w(I_Pracw  ) & ! [prod] c->r
           - w(I_Prevp  ) & ! [loss] r->v
           + ( &
           + w(I_Psacw  ) & ! [prod] c->s
           + w(I_Pgacw  ) & ! [prod] c->g
           + w(I_Psmlt  ) & ! [prod] s->r
           + w(I_Pgmlt  ) & ! [prod] g->r
           ) * ( 1.0_RP - ice ) &
           + ( &
           - w(I_Piacr_s) & ! [loss] r->s
           - w(I_Psacr_s) & ! [loss] r->s
           - w(I_Piacr_g) & ! [loss] r->g
           - w(I_Psacr_g) & ! [loss] r->g
           - w(I_Pgacr  ) & ! [loss] r->g
           - w(I_Pgfrz  ) & ! [loss] r->g
           ) * ice

       fac_sw = 0.5_RP + sign( 0.5_RP, net+EPS ) ! if production > loss , fac_sw=1
       fac    = (          fac_sw ) &
              + ( 1.0_RP - fac_sw ) * min( -w(I_dqr_dt)/(net-fac_sw), 1.0_RP ) ! loss limiter

       fac_warm = fac * (1.0_RP-ice) + 1.0_RP * ice
       fac_ice  = fac * ice + 1.0_RP * (1.0_RP-ice)

       w(I_Praut) = w(I_Praut) * fac
       w(I_Pracw) = w(I_Pracw) * fac
       w(I_Prevp) = w(I_Prevp) * fac

       w(I_Psacw) = w(I_Psacw) * fac_warm
       w(I_Pgacw) = w(I_Pgacw) * fac_warm
       w(I_Psmlt) = w(I_Psmlt) * fac_warm
       w(I_Pgmlt) = w(I_Pgmlt) * fac_warm

       w(I_Piacr_s) = w(I_Piacr_s) * fac_ice
       w(I_Psacr_s) = w(I_Psacr_s) * fac_ice
       w(I_Piacr_g) = w(I_Piacr_g) * fac_ice
       w(I_Psacr_g) = w(I_Psacr_g) * fac_ice
       w(I_Pgacr  ) = w(I_Pgacr  ) * fac_ice
       w(I_Pgfrz  ) = w(I_Pgfrz  ) * fac_ice

       !### in the case of Sice-1 > 0, limiter is applied to qv before qs,qg
       ! [QV]
       net = w(I_Prevp) & ! [prod] r->v
           + w(I_Pssub) & ! [prod] s->v
           + w(I_Pgsub) & ! [prod] g->v
           - w(I_Psdep) & ! [loss] v->s
           - w(I_Pgdep)   ! [loss] v->g

       fac_sw = 0.5_RP + sign( 0.5_RP, net+EPS ) ! if production > loss , fac_sw=1
       fac    = (          fac_sw ) &
              + ( 1.0_RP - fac_sw ) * min( -w(I_dqv_dt)/(net-fac_sw), 1.0_RP ) ! loss limiter

       fac    = (          w(I_delta3) ) * fac &
              + ( 1.0_RP - w(I_delta3) )

       fac = fac * ice + 1.0_RP * (1.0_RP - ice)

       w(I_Prevp) = w(I_Prevp) * fac
       w(I_Pssub) = w(I_Pssub) * fac
       w(I_Pgsub) = w(I_Pgsub) * fac
       w(I_Psdep) = w(I_Psdep) * fac
       w(I_Pgdep) = w(I_Pgdep) * fac

       ! [QS]
       net = &
           - w(I_Pgacs  ) & ! [loss] s->g
           + ( &
           - w(I_Psmlt  ) & ! [loss] s->r
           ) * ( 1.0_RP - ice ) &
           + ( &
           + w(I_Psdep  ) & ! [prod] v->s
           + w(I_Psacw  ) & ! [prod] c->s
           + w(I_Psfw   ) & ! [prod] c->s
           + w(I_Piacr_s) & ! [prod] r->s
           + w(I_Psacr_s) & ! [prod] r->s
           + w(I_Psaut  ) & ! [prod] i->s
           + w(I_Praci_s) & ! [prod] i->s
           + w(I_Psaci  ) & ! [prod] i->s
           + w(I_Psfi   ) & ! [prod] i->s
           - w(I_Pssub  ) & ! [loss] s->v
           - w(I_Pgaut  ) & ! [loss] s->g
           - w(I_Pracs  ) & ! [loss] s->g
           )

       fac_sw = 0.5_RP + sign( 0.5_RP, net+EPS ) ! if production > loss , fac_sw=1
       fac    = (          fac_sw ) &
              + ( 1.0_RP - fac_sw ) * min( -w(I_dqs_dt)/(net-fac_sw), 1.0_RP ) ! loss limiter

       fac_warm = fac * (1.0_RP - ice) + 1.0_RP * ice
       fac_ice  = fac * ice + 1.0_RP * (1.0_RP-ice)

       w(I_Pgacs) = w(I_Pgacs) * fac

       w(I_Psmlt) = w(I_Psmlt) * fac_warm

       w(I_Psdep  ) = w(I_Psdep  ) * fac_ice
       w(I_Psacw  ) = w(I_Psacw  ) * fac_ice
       w(I_Psfw   ) = w(I_Psfw   ) * fac_ice
       w(I_Piacr_s) = w(I_Piacr_s) * fac_ice
       w(I_Psacr_s) = w(I_Psacr_s) * fac_ice
       w(I_Psaut  ) = w(I_Psaut  ) * fac_ice
       w(I_Praci_s) = w(I_Praci_s) * fac_ice
       w(I_Psaci  ) = w(I_Psaci  ) * fac_ice
       w(I_Psfi   ) = w(I_Psfi   ) * fac_ice
       w(I_Pssub  ) = w(I_Pssub  ) * fac_ice
       w(I_Pracs  ) = w(I_Pracs  ) * fac_ice
       w(I_Pgaut  ) = w(I_Pgaut  ) * fac_ice

       ! [QG]
       net = w(I_Pgacs  ) & ! [prod] s->g
           + ( &
           - w(I_Pgmlt  ) & ! [loss] g->r
           ) * ( 1.0_RP - ice ) &
           + ( &
           + w(I_Pgdep  ) & ! [prod] v->g
           + w(I_Pgacw  ) & ! [prod] c->g
           + w(I_Piacr_g) & ! [prod] r->g
           + w(I_Psacr_g) & ! [prod] r->g
           + w(I_Pgacr  ) & ! [prod] r->g
           + w(I_Pgfrz  ) & ! [prod] r->g
           + w(I_Praci_g) & ! [prod] i->g
           + w(I_Pgaci  ) & ! [prod] i->g
           + w(I_Pgaut  ) & ! [prod] s->g
           + w(I_Pracs  ) & ! [prod] s->g
           - w(I_Pgsub  ) & ! [loss] g->v
           ) * ice

       fac_sw = 0.5_RP + sign( 0.5_RP, net+EPS ) ! if production > loss , fac_sw=1
       fac    = (          fac_sw ) &
              + ( 1.0_RP - fac_sw ) * min( -w(I_dqg_dt)/(net-fac_sw), 1.0_RP ) ! loss limiter

       fac_warm = fac * (1.0_RP - ice) + 1.0_RP * ice
       fac_ice  = fac * ice + 1.0_RP * (1.0_RP-ice)

       w(I_Pgacs) = w(I_Pgacs) * fac

       w(I_Pgmlt) = w(I_Pgmlt) * fac_warm

       w(I_Pgdep  ) = w(I_Pgdep  ) * fac_ice
       w(I_Pgacw  ) = w(I_Pgacw  ) * fac_ice
       w(I_Piacr_g) = w(I_Piacr_g) * fac_ice
       w(I_Psacr_g) = w(I_Psacr_g) * fac_ice
       w(I_Pgacr  ) = w(I_Pgacr  ) * fac_ice
       w(I_Pgfrz  ) = w(I_Pgfrz  ) * fac_ice
       w(I_Praci_g) = w(I_Praci_g) * fac_ice
       w(I_Pgaci  ) = w(I_Pgaci  ) * fac_ice
       w(I_Pgaut  ) = w(I_Pgaut  ) * fac_ice
       w(I_Pracs  ) = w(I_Pracs  ) * fac_ice
       w(I_Pgsub  ) = w(I_Pgsub  ) * fac_ice



       ! re-calc net production & loss
       tend(I_QC) = &
                  - w(I_Praut  ) & ! [loss] c->r
                  - w(I_Pracw  ) & ! [loss] c->r
                  - w(I_Psacw  ) & ! [loss] c->s
                  - w(I_Pgacw  ) & ! [loss] c->g
                  - w(I_Psfw   ) * ice

       tend(I_QR) = w(I_Praut  ) & ! [prod] c->r
                  + w(I_Pracw  ) & ! [prod] c->r
                  - w(I_Prevp  ) & ! [loss] r->v
                  + ( &
                  + w(I_Psacw  ) & ! [prod] c->s=r
                  + w(I_Pgacw  ) & ! [prod] c->g=r
                  + w(I_Psmlt  ) & ! [prod] s->r
                  + w(I_Pgmlt  ) & ! [prod] g->r
                  ) * ( 1.0_RP - ice ) &
                  + ( &
                  - w(I_Piacr_s) & ! [loss] r->s
                  - w(I_Psacr_s) & ! [loss] r->s
                  - w(I_Piacr_g) & ! [loss] r->g
                  - w(I_Psacr_g) & ! [loss] r->g
                  - w(I_Pgacr  ) & ! [loss] r->g
                  - w(I_Pgfrz  ) & ! [loss] r->g
                  ) * ice

       tend(I_QI) = ( &
                  - w(I_Psaut  ) & ! [loss] i->s
                  - w(I_Praci_s) & ! [loss] i->s
                  - w(I_Psaci  ) & ! [loss] i->s
                  - w(I_Psfi   ) & ! [loss] i->s
                  - w(I_Praci_g) & ! [loss] i->g
                  - w(I_Pgaci  ) & ! [loss] i->g
                  ) * ice

       tend(I_QS) = &
                  - w(I_Pgacs  ) & ! [loss] s->g
                  + ( &
                  - w(I_Psmlt  ) & ! [loss] s->r
                  ) * ( 1.0_RP - ice ) &
                  + ( &
                  + w(I_Psdep  ) & ! [prod] v->s
                  + w(I_Psacw  ) & ! [prod] c->s
                  + w(I_Psfw   ) & ! [prod] c->s
                  + w(I_Piacr_s) & ! [prod] r->s
                  + w(I_Psacr_s) & ! [prod] r->s
                  + w(I_Psaut  ) & ! [prod] i->s
                  + w(I_Psaci  ) & ! [prod] i->s
                  + w(I_Praci_s) & ! [prod] i->s
                  + w(I_Psfi   ) & ! [prod] i->s
                  - w(I_Pssub  ) & ! [loss] s->v
                  - w(I_Pgaut  ) & ! [loss] s->g
                  - w(I_Pracs  ) & ! [loss] s->g
                  ) * ice

       tend(I_QG) = w(I_Pgacs  ) & ! [prod] s->g
                  + ( &
                  - w(I_Pgmlt  ) & ! [loss] g->r
                  ) * ( 1.0_RP - ice ) &
                  + ( &
                  + w(I_Pgdep  ) & ! [prod] v->g
                  + w(I_Pgacw  ) & ! [prod] c->g
                  + w(I_Piacr_g) & ! [prod] r->g
                  + w(I_Psacr_g) & ! [prod] r->g
                  + w(I_Pgacr  ) & ! [prod] r->g
                  + w(I_Pgfrz  ) & ! [prod] r->g
                  + w(I_Praci_g) & ! [prod] i->g
                  + w(I_Pgaci  ) & ! [prod] i->g
                  + w(I_Pgaut  ) & ! [prod] s->g
                  + w(I_Pracs  ) & ! [prod] s->g
                  - w(I_Pgsub  ) & ! [loss] g->v
                  ) * ice

       tend(I_QC) = max( tend(I_QC), -w(I_dqc_dt) )
       tend(I_QR) = max( tend(I_QR), -w(I_dqr_dt) )
       tend(I_QI) = max( tend(I_QI), -w(I_dqi_dt) )
       tend(I_QS) = max( tend(I_QS), -w(I_dqs_dt) )
       tend(I_QG) = max( tend(I_QG), -w(I_dqg_dt) )

       tend(I_QV) = - ( tend(I_QC) &
                      + tend(I_QR) &
                      + tend(I_QI) &
                      + tend(I_QS) &
                      + tend(I_QG) )

       QTRC_t(k,i,j,I_QV) = tend(I_QV)
       QTRC_t(k,i,j,I_QC) = tend(I_QC)
       QTRC_t(k,i,j,I_QR) = tend(I_QR)
       QTRC_t(k,i,j,I_QI) = tend(I_QI)
       QTRC_t(k,i,j,I_QS) = tend(I_QS)
       QTRC_t(k,i,j,I_QG) = tend(I_QG)

       QTRC0(k,i,j,I_QV) = QTRC0(k,i,j,I_QV) + QTRC_t(k,i,j,I_QV) * dt
       QTRC0(k,i,j,I_QC) = QTRC0(k,i,j,I_QC) + QTRC_t(k,i,j,I_QC) * dt
       QTRC0(k,i,j,I_QR) = QTRC0(k,i,j,I_QR) + QTRC_t(k,i,j,I_QR) * dt
       QTRC0(k,i,j,I_QI) = QTRC0(k,i,j,I_QI) + QTRC_t(k,i,j,I_QI) * dt
       QTRC0(k,i,j,I_QS) = QTRC0(k,i,j,I_QS) + QTRC_t(k,i,j,I_QS) * dt
       QTRC0(k,i,j,I_QG) = QTRC0(k,i,j,I_QG) + QTRC_t(k,i,j,I_QG) * dt

!      w3d(k,i,j,ip) = w(ip)
    enddo
    enddo
    enddo

    !-----< ice generation process >-----
    if ( MP_doexpricit_icegen ) then
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          ! store to work
          dens = DENS0(k,i,j)
          temp = TEMP0(k,i,j)
          pres = PRES0(k,i,j)

          ! saturation ratio S
          Sice = QTRC0(k,i,j,I_QV) / max( QSATI(k,i,j), EPS )

          Rdens = 1.0_RP / dens
          temc  = temp - TEM00
          rhoqi = dens * QTRC0(k,i,j,I_QI)

          Da = ( Da0 + dDa_dT * temc )
          Kd = ( Dw0 + dDw_dT * temc ) * PRE00 / pres

          Giv  = 1.0_RP / ( LHS0/(Da*temp) * ( LHS0/(Rvap*temp) - 1.0_RP ) + 1.0_RP/(Kd*dens*QSATI(k,i,j)) )

          ! [Pidep] deposition/sublimation : v->i or i->v
          sw = ( 0.5_RP + sign(0.5_RP, 0.0_RP - temc ) ) &
             * ( 0.5_RP + sign(0.5_RP, rhoqi         ) ) ! if T < 0C & ice exists, sw=1

          rhoqi = max(rhoqi,EPS)

          XNi = min( max( 5.38E7_RP * rhoqi**0.75_RP, 1.E3_RP ), 1.E6_RP )
          XMi = rhoqi / XNi
          Di  = min( Di_a * sqrt(XMi), Di_max )

          Pidep = sw * 4.0_RP * Di * XNi * Rdens * ( Sice-1.0_RP ) * Giv

          ! [Pigen] ice nucleation : v->i
          sw = ( 0.5_RP + sign(0.5_RP, 0.0_RP - temc ) ) &
             * ( 0.5_RP + sign(0.5_RP, Sice - 1.0_RP ) ) ! if T < 0C & Satulated, sw=1

          Ni0 = max( exp(-0.1_RP*temc), 1.0_RP ) * 1000.0_RP
          Qi0 = 4.92E-11 * Ni0**1.33_RP * Rdens

          Pigen = sw * max( min( Qi0-QTRC0(k,i,j,I_QI), QTRC0(k,i,j,I_QV)-QSATI(k,i,j) ), 0.0_RP ) / dt

          ! update Qv, Qi with limiter
          dq   = ( Pigen + Pidep ) * dt
          qtmp = QTRC0(k,i,j,I_QV) - dq
          if ( dq > 0.0_RP ) then ! v->i
             qtmp = max( qtmp, QSATI(k,i,j) )
          else                    ! i->v
             qtmp = min( qtmp, QSATI(k,i,j) )
          endif
          dq   = QTRC0(k,i,j,I_QV) - qtmp
          qtmp = QTRC0(k,i,j,I_QI) + dq
          qtmp = max( qtmp, 0.0_RP )
          dq   = qtmp - QTRC0(k,i,j,I_QI)
          QTRC0 (k,i,j,I_QI) = QTRC0 (k,i,j,I_QI) + dq
          QTRC0 (k,i,j,I_QV) = QTRC0 (k,i,j,I_QV) - dq
          QTRC_t(k,i,j,I_QI) = QTRC_t(k,i,j,I_QI) + dq /dt
          QTRC_t(k,i,j,I_QV) = QTRC_t(k,i,j,I_QV) - dq /dt

          ! [Pifrz/Pimlt] freezing and melting

          ! homogenious freezing : c->i
          sw = ( 0.5_RP - sign(0.5_RP, temc + 40.0_RP ) ) ! if T < -40C, sw=1

          dq = sw * QTRC0(k,i,j,I_qc)
          QTRC0 (k,i,j,I_QI) = QTRC0 (k,i,j,I_QI) + dq
          QTRC0 (k,i,j,I_QC) = QTRC0 (k,i,j,I_QC) - dq
          QTRC_t(k,i,j,I_QI) = QTRC_t(k,i,j,I_QI) + dq /dt
          QTRC_t(k,i,j,I_QC) = QTRC_t(k,i,j,I_QC) - dq /dt

          ! heteroginous freezing : c->i
          sw = ( 0.5_RP + sign(0.5_RP, temc + 40.0_RP ) ) &
             * ( 0.5_RP + sign(0.5_RP, 0.0_RP -  temc ) ) ! if -40C < T < 0C, sw=1

          dq = sw * ( dens / DWATR * QTRC0(k,i,j,I_QC)**2 / (Nc(k,i,j)*1.E+6) ) &
                  * B_frz * ( exp(-A_frz*temc) - 1.0_RP ) * dt
          qtmp = QTRC0(k,i,j,I_QC) - dq
          qtmp = max( qtmp, 0.0_RP )
          dq   = QTRC0(k,i,j,I_QC) - qtmp
          QTRC0 (k,i,j,I_QI) = QTRC0 (k,i,j,I_QI) + dq
          QTRC0 (k,i,j,I_QC) = QTRC0 (k,i,j,I_QC) - dq
          QTRC_t(k,i,j,I_QI) = QTRC_t(k,i,j,I_QI) + dq /dt
          QTRC_t(k,i,j,I_QC) = QTRC_t(k,i,j,I_QC) - dq /dt

          ! melting : i->c
          sw = ( 0.5_RP - sign(0.5_RP, 0.0_RP -  temc ) ) ! if 0C < T, sw=1

          dq = - sw * QTRC0(k,i,j,I_QI)
          QTRC0 (k,i,j,I_QI) = QTRC0 (k,i,j,I_QI) + dq
          QTRC0 (k,i,j,I_QC) = QTRC0 (k,i,j,I_QC) - dq
          QTRC_t(k,i,j,I_QI) = QTRC_t(k,i,j,I_QI) + dq /dt
          QTRC_t(k,i,j,I_QC) = QTRC_t(k,i,j,I_QC) - dq /dt

       enddo
       enddo
       enddo
    endif

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       RHOE_t(k,i,j) = - DENS0(k,i,j) * ( LHV * QTRC_t(k,i,j,I_QV) &
                                        - LHF * QTRC_t(k,i,j,I_QI) &
                                        - LHF * QTRC_t(k,i,j,I_QS) &
                                        - LHF * QTRC_t(k,i,j,I_QG) )

       RHOE0(k,i,j) = RHOE0(k,i,j) + RHOE_t(k,i,j) * dt
    enddo
    enddo
    enddo

!    if( IO_L ) write(IO_FID_LOG,*) "QV      MAX/MIN:", maxval(QTRC0 (:,:,:,I_QV)), minval(QTRC0 (:,:,:,I_QV))
!    if( IO_L ) write(IO_FID_LOG,*) "QC      MAX/MIN:", maxval(QTRC0 (:,:,:,I_QC)), minval(QTRC0 (:,:,:,I_QC))
!    if( IO_L ) write(IO_FID_LOG,*) "QR      MAX/MIN:", maxval(QTRC0 (:,:,:,I_QR)), minval(QTRC0 (:,:,:,I_QR))
!    if( IO_L ) write(IO_FID_LOG,*) "QI      MAX/MIN:", maxval(QTRC0 (:,:,:,I_QI)), minval(QTRC0 (:,:,:,I_QI))
!    if( IO_L ) write(IO_FID_LOG,*) "QS      MAX/MIN:", maxval(QTRC0 (:,:,:,I_QS)), minval(QTRC0 (:,:,:,I_QS))
!    if( IO_L ) write(IO_FID_LOG,*) "QG      MAX/MIN:", maxval(QTRC0 (:,:,:,I_QG)), minval(QTRC0 (:,:,:,I_QG))
!    if( IO_L ) write(IO_FID_LOG,*) "QV_t    MAX/MIN:", maxval(QTRC_t(:,:,:,I_QV)), minval(QTRC_t(:,:,:,I_QV))
!    if( IO_L ) write(IO_FID_LOG,*) "QC_t    MAX/MIN:", maxval(QTRC_t(:,:,:,I_QC)), minval(QTRC_t(:,:,:,I_QC))
!    if( IO_L ) write(IO_FID_LOG,*) "QR_t    MAX/MIN:", maxval(QTRC_t(:,:,:,I_QR)), minval(QTRC_t(:,:,:,I_QR))
!    if( IO_L ) write(IO_FID_LOG,*) "QI_t    MAX/MIN:", maxval(QTRC_t(:,:,:,I_QI)), minval(QTRC_t(:,:,:,I_QI))
!    if( IO_L ) write(IO_FID_LOG,*) "QS_t    MAX/MIN:", maxval(QTRC_t(:,:,:,I_QS)), minval(QTRC_t(:,:,:,I_QS))
!    if( IO_L ) write(IO_FID_LOG,*) "QG_t    MAX/MIN:", maxval(QTRC_t(:,:,:,I_QG)), minval(QTRC_t(:,:,:,I_QG))
!    if( IO_L ) write(IO_FID_LOG,*) "RHOE0   MAX/MIN:", maxval(RHOE0 (:,:,:)),      minval(RHOE0 (:,:,:))
!    if( IO_L ) write(IO_FID_LOG,*) "RHOE_t  MAX/MIN:", maxval(RHOE_t(:,:,:)),      minval(RHOE_t(:,:,:))
!    do ip = 1, w_nmax
!       if( IO_L ) write(IO_FID_LOG,*) w_name(ip), "MAX/MIN:", &
!                  maxval(w3d(KS:KE,IS:IE,JS:JE,ip)), minval(w3d(KS:KE,IS:IE,JS:JE,ip))
!
!       call HIST_in( w_name(ip), w3d(:,:,:,ip), 'individual tendency term in tomita08', 'kg/kg/s' )
!    enddo

    call PROF_rapend  ('MP_tomita08', 3)

    return
  end subroutine MP_tomita08

  !-----------------------------------------------------------------------------
  !> Lin-type cold rain microphysics (terminal velocity)
  subroutine MP_tomita08_vterm( &
       vterm, &
       DENS0, &
       TEMP0, &
       QTRC0  )
    use scale_const, only: &
       TEM00 => CONST_TEM00
    use scale_tracer, only: &
       QA
    implicit none

    real(RP), intent(out) :: vterm(KA,IA,JA,QA_MP-1)
    real(RP), intent(in)  :: DENS0(KA,IA,JA)
    real(RP), intent(in)  :: TEMP0(KA,IA,JA)
    real(RP), intent(in)  :: QTRC0(KA,IA,JA,QA)

    real(RP) :: dens
    real(RP) :: temc
    real(RP) :: q(QA)

    real(RP) :: rho_fact ! density factor

    real(RP) :: RLMDr, RLMDs, RLMDg
    real(RP) :: RLMDr_dr, RLMDs_ds, RLMDg_dg

    !---< Roh and Satoh (2014) >---
    real(RP) :: tems, rhoqs, Xs2
    real(RP) :: MOMs_0bs, RMOMs_Vt
    real(RP) :: coef_at(4), coef_bt(4)
    real(RP) :: loga_, b_, nm

    real(RP) :: zerosw
    integer  :: k, i, j, iq
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       ! store to work
       dens = DENS0(k,i,j)
       temc = TEMP0(k,i,j) - TEM00
       do iq = QS_MP, QE_MP
          q(iq) = QTRC0(k,i,j,iq)
       enddo

       rho_fact = sqrt( dens00 / dens )

       ! slope parameter lambda (Rain)
       zerosw = 0.5_RP - sign(0.5_RP, q(I_QR) - 1.E-12_RP )
       RLMDr  = sqrt(sqrt( dens * q(I_QR) / ( Ar * N0r * GAM_1br ) + zerosw )) * ( 1.0_RP - zerosw )

       RLMDr_dr = sqrt( RLMDr )       ! **Dr

       ! slope parameter lambda (Snow)
       zerosw = 0.5_RP - sign(0.5_RP, q(I_QS) - 1.E-12_RP )
       RLMDs  = sqrt(sqrt( dens * q(I_QS) / ( As * N0s * GAM_1bs ) + zerosw )) * ( 1.0_RP - zerosw )

       RLMDs_ds = sqrt( sqrt(RLMDs) ) ! **Ds
       RMOMs_Vt = GAM_1bsds / GAM_1bs * RLMDs_ds

       !---< modification by Roh and Satoh (2014) >---
       ! bimodal size distribution of snow
       rhoqs  = dens * q(I_QS)
       zerosw = 0.5_RP - sign(0.5_RP, rhoqs - 1.E-12_RP )
       Xs2    = rhoqs / As * ( 1.0_RP - zerosw )

       tems       = min( -0.1_RP, temc )
       coef_at(1) = coef_a( 1) + tems * ( coef_a( 2) + tems * ( coef_a( 5) + tems * coef_a( 9) ) )
       coef_at(2) = coef_a( 3) + tems * ( coef_a( 4) + tems *   coef_a( 7) )
       coef_at(3) = coef_a( 6) + tems *   coef_a( 8)
       coef_at(4) = coef_a(10)
       coef_bt(1) = coef_b( 1) + tems * ( coef_b( 2) + tems * ( coef_b( 5) + tems * coef_b( 9) ) )
       coef_bt(2) = coef_b( 3) + tems * ( coef_b( 4) + tems *   coef_b( 7) )
       coef_bt(3) = coef_b( 6) + tems *   coef_b( 8)
       coef_bt(4) = coef_b(10)
       ! 0 + Bs(=2) moment
       nm = 2.0_RP
       loga_  = coef_at(1) + nm * ( coef_at(2) + nm * ( coef_at(3) + nm * coef_at(4) ) )
          b_  = coef_bt(1) + nm * ( coef_bt(2) + nm * ( coef_bt(3) + nm * coef_bt(4) ) )
       MOMs_0bs = 10.0_RP**loga_ * Xs2**b_
       ! Bs(=2) + Ds(=0.25) moment
       nm = 2.25_RP
       loga_  = coef_at(1) + nm * ( coef_at(2) + nm * ( coef_at(3) + nm * coef_at(4) ) )
          b_  = coef_bt(1) + nm * ( coef_bt(2) + nm * ( coef_bt(3) + nm * coef_bt(4) ) )
       RMOMs_Vt = (        sw_roh2014 ) * 10.0_RP**loga_ * Xs2**b_ / ( MOMs_0bs + zerosw ) &
                + ( 1.0_RP-sw_roh2014 ) * RMOMs_Vt

       ! slope parameter lambda (Graupel)
       zerosw = 0.5_RP - sign(0.5_RP, q(I_QG) - 1.E-12_RP )
       RLMDg  = sqrt(sqrt( dens * q(I_QG) / ( Ag * N0g * GAM_1bg ) + zerosw )) * ( 1.0_RP - zerosw )

       RLMDg_dg = sqrt( RLMDg )       ! **Dg

       !---< terminal velocity >
       zerosw = 0.5_RP + sign(0.5_RP, q(I_QI) - 1.E-8_RP )
       vterm(k,i,j,I_mp_QC) = 0.0_RP
       vterm(k,i,j,I_mp_QI) = -3.29_RP * ( dens * q(I_QI) * zerosw )**0.16_RP
       vterm(k,i,j,I_mp_QR) = -Cr * rho_fact * GAM_1brdr / GAM_1br * RLMDr_dr
       vterm(k,i,j,I_mp_QS) = -Cs * rho_fact * RMOMs_Vt
       vterm(k,i,j,I_mp_QG) = -Cg * rho_fact * GAM_1bgdg / GAM_1bg * RLMDg_dg
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       vterm(    1:KS-1,i,j,:) = 0.0_RP
       vterm(KE+1:KA   ,i,j,:) = 0.0_RP
    enddo
    enddo

    return
  end subroutine MP_tomita08_vterm

  !-----------------------------------------------------------------------------
  subroutine MP_tomita08_BergeronParam( &
       temp, &
       a1,   &
       a2,   &
       ma2   )
    use scale_const, only: &
       TEM00 => CONST_TEM00
    implicit none

    real(RP), intent(in)  :: temp(KA,IA,JA)
    real(RP), intent(out) :: a1  (KA,IA,JA)
    real(RP), intent(out) :: a2  (KA,IA,JA)
    real(RP), intent(out) :: ma2 (KA,IA,JA)

    real(RP) :: a1_tab(32)
    real(RP) :: a2_tab(32)

    data a1_tab / 0.0001E-7_RP, 0.7939E-7_RP, 0.7841E-6_RP, 0.3369E-5_RP, 0.4336E-5_RP, &
                  0.5285E-5_RP, 0.3728E-5_RP, 0.1852E-5_RP, 0.2991E-6_RP, 0.4248E-6_RP, &
                  0.7434E-6_RP, 0.1812E-5_RP, 0.4394E-5_RP, 0.9145E-5_RP, 0.1725E-4_RP, &
                  0.3348E-4_RP, 0.1725E-4_RP, 0.9175E-5_RP, 0.4412E-5_RP, 0.2252E-5_RP, &
                  0.9115E-6_RP, 0.4876E-6_RP, 0.3473E-6_RP, 0.4758E-6_RP, 0.6306E-6_RP, &
                  0.8573E-6_RP, 0.7868E-6_RP, 0.7192E-6_RP, 0.6513E-6_RP, 0.5956E-6_RP, &
                  0.5333E-6_RP, 0.4834E-6_RP  /

    data a2_tab / 0.0100_RP, 0.4006_RP, 0.4831_RP, 0.5320_RP, 0.5307_RP, &
                  0.5319_RP, 0.5249_RP, 0.4888_RP, 0.3849_RP, 0.4047_RP, &
                  0.4318_RP, 0.4771_RP, 0.5183_RP, 0.5463_RP, 0.5651_RP, &
                  0.5813_RP, 0.5655_RP, 0.5478_RP, 0.5203_RP, 0.4906_RP, &
                  0.4447_RP, 0.4126_RP, 0.3960_RP, 0.4149_RP, 0.4320_RP, &
                  0.4506_RP, 0.4483_RP, 0.4460_RP, 0.4433_RP, 0.4413_RP, &
                  0.4382_RP, 0.4361_RP  /

    real(RP) :: temc
    integer  :: itemc
    real(RP) :: fact

    integer :: k, i, j
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       temc  = min( max( temp(k,i,j)-TEM00, -30.99_RP ), 0.0_RP ) ! 0C <= T  <  31C
       itemc = int( -temc ) + 1                                   ! 1  <= iT <= 31
       fact  = - ( temc + real(itemc-1,kind=8) )

       a1(k,i,j) = ( 1.0_RP-fact ) * a1_tab(itemc  ) &
                 + (        fact ) * a1_tab(itemc+1)

       a2(k,i,j) = ( 1.0_RP-fact ) * a2_tab(itemc  ) &
                 + (        fact ) * a2_tab(itemc+1)

       ma2(k,i,j) = 1.0_RP - a2(k,i,j)

       a1(k,i,j) = a1(k,i,j) * 1.E-3_RP**ma2(k,i,j) ! [g->kg]
    enddo
    enddo
    enddo

    return
  end subroutine MP_tomita08_BergeronParam

  !-----------------------------------------------------------------------------
  !> Calculate Cloud Fraction
  subroutine ATMOS_PHY_MP_tomita08_CloudFraction( &
       cldfrac, &
       QTRC     )
    use scale_grid_index
    use scale_tracer, only: &
       QA
    implicit none

    real(RP), intent(out) :: cldfrac(KA,IA,JA)
    real(RP), intent(in)  :: QTRC   (KA,IA,JA,QA)

    real(RP) :: qcriteria = 0.005E-3_RP ! 0.005g/kg, Tompkins & Craig

    real(RP) :: qhydro
    integer  :: k, i, j, iq
    !---------------------------------------------------------------------------

    do j  = JS, JE
    do i  = IS, IE
    do k  = KS, KE
       qhydro = 0.0_RP
       do iq = QS_MP+1, QE_MP
          qhydro = qhydro + QTRC(k,i,j,iq)
       enddo
       cldfrac(k,i,j) = 0.5_RP + sign(0.5_RP,qhydro-qcriteria)
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_MP_tomita08_CloudFraction

  !-----------------------------------------------------------------------------
  !> Calculate Effective Radius
  subroutine ATMOS_PHY_MP_tomita08_EffectiveRadius( &
       Re,    &
       QTRC0, &
       DENS0, &
       TEMP0  )
    use scale_grid_index
    use scale_const, only: &
       TEM00 => CONST_TEM00
    use scale_tracer, only: &
       QA
    use scale_atmos_hydrometeor, only: &
       N_HYD
    implicit none

    real(RP), intent(out) :: Re   (KA,IA,JA,N_HYD) ! effective radius          [cm]
    real(RP), intent(in)  :: QTRC0(KA,IA,JA,QA)    ! tracer mass concentration [kg/kg]
    real(RP), intent(in)  :: DENS0(KA,IA,JA)       ! density                   [kg/m3]
    real(RP), intent(in)  :: TEMP0(KA,IA,JA)       ! temperature               [K]

    real(RP) :: dens
    real(RP) :: temc
    real(RP) :: RLMDr, RLMDs, RLMDg

    real(RP), parameter :: um2cm = 100.0_RP

    !---< Roh and Satoh (2014) >---
    real(RP) :: tems, rhoqs, Xs2
    real(RP) :: MOMs_2, MOMs_1bs
    real(RP) :: coef_at(4), coef_bt(4)
    real(RP) :: loga_, b_, nm

    real(RP) :: zerosw
    integer  :: k, i, j, iq
    !---------------------------------------------------------------------------

    Re(:,:,:,I_HC) =  re_qc * um2cm
    Re(:,:,:,I_HI) =  re_qi * um2cm
    Re(:,:,:,I_HG+1:) = 0.0_RP

    do j  = JS, JE
    do i  = IS, IE
    do k  = KS, KE
       dens = DENS0(k,i,j)
       temc = TEMP0(k,i,j) - TEM00

       ! slope parameter lambda
       zerosw = 0.5_RP - sign(0.5_RP, QTRC0(k,i,j,I_QR) - 1.E-12_RP )
       RLMDr = sqrt(sqrt( dens * QTRC0(k,i,j,I_QR) / ( Ar * N0r * GAM_1br ) + zerosw )) * ( 1.0_RP - zerosw )
       ! Effective radius is defined by r3m/r2m = 1.5/lambda
       Re(k,i,j,I_HR) = 1.5_RP * RLMDr * um2cm


       zerosw = 0.5_RP - sign(0.5_RP, QTRC0(k,i,j,I_QS) - 1.E-12_RP )
       RLMDs = sqrt(sqrt( dens * QTRC0(k,i,j,I_QS) / ( As * N0s * GAM_1bs ) + zerosw )) * ( 1.0_RP - zerosw )
       !---< modification by Roh and Satoh (2014) >---
       ! bimodal size distribution of snow
       rhoqs  = dens * QTRC0(k,i,j,I_QS)
       zerosw = 0.5_RP - sign(0.5_RP, rhoqs - 1.E-12_RP )
       Xs2    = rhoqs / As * ( 1.0_RP - zerosw )
       tems       = min( -0.1_RP, temc )
       coef_at(1) = coef_a( 1) + tems * ( coef_a( 2) + tems * ( coef_a( 5) + tems * coef_a( 9) ) )
       coef_at(2) = coef_a( 3) + tems * ( coef_a( 4) + tems *   coef_a( 7) )
       coef_at(3) = coef_a( 6) + tems *   coef_a( 8)
       coef_at(4) = coef_a(10)
       coef_bt(1) = coef_b( 1) + tems * ( coef_b( 2) + tems * ( coef_b( 5) + tems * coef_b( 9) ) )
       coef_bt(2) = coef_b( 3) + tems * ( coef_b( 4) + tems *   coef_b( 7) )
       coef_bt(3) = coef_b( 6) + tems *   coef_b( 8)
       coef_bt(4) = coef_b(10)
       ! 2nd moment
       MOMs_2 = Xs2
       ! 1 + Bs(=2) moment
       nm = 3.0_RP
       loga_  = coef_at(1) + nm * ( coef_at(2) + nm * ( coef_at(3) + nm * coef_at(4) ) )
          b_  = coef_bt(1) + nm * ( coef_bt(2) + nm * ( coef_bt(3) + nm * coef_bt(4) ) )
       MOMs_1bs = 10.0_RP**loga_ * Xs2**b_
       Re(k,i,j,I_HS) = (        sw_roh2014 ) * 0.5_RP * MOMs_1bs / ( MOMs_2 + zerosw ) * um2cm &
                      + ( 1.0_RP-sw_roh2014 ) * 1.5_RP * RLMDs * um2cm

       zerosw = 0.5_RP - sign(0.5_RP, QTRC0(k,i,j,I_QG) - 1.E-12_RP )
       RLMDg = sqrt(sqrt( dens * QTRC0(k,i,j,I_QG) / ( Ag * N0g * GAM_1bg ) + zerosw )) * ( 1.0_RP - zerosw )
       Re(k,i,j,I_HG) = 1.5_RP * RLMDg * um2cm
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_MP_tomita08_EffectiveRadius

  !-----------------------------------------------------------------------------
  !> Calculate mixing ratio of each category
  subroutine ATMOS_PHY_MP_tomita08_Mixingratio( &
       Qe,   &
       QTRC0 )
    use scale_grid_index
    use scale_tracer, only: &
       QA
    use scale_atmos_hydrometeor, only: &
       N_HYD
    implicit none

    real(RP), intent(out) :: Qe   (KA,IA,JA,N_HYD) ! mixing ratio of each cateory [kg/kg]
    real(RP), intent(in)  :: QTRC0(KA,IA,JA,QA)    ! tracer mass concentration [kg/kg]

    integer :: ihydro, iqa
    !---------------------------------------------------------------------------

    Qe(:,:,:,I_HC) = QTRC0(:,:,:,I_QC)
    Qe(:,:,:,I_HR) = QTRC0(:,:,:,I_QR)
    Qe(:,:,:,I_HI) = QTRC0(:,:,:,I_QI)
    Qe(:,:,:,I_HS) = QTRC0(:,:,:,I_QS)
    Qe(:,:,:,I_HG) = QTRC0(:,:,:,I_QG)
    Qe(:,:,:,I_HG+1:) = 0.0_RP

    return
  end subroutine ATMOS_PHY_MP_tomita08_Mixingratio

end module scale_atmos_phy_mp_tomita08
