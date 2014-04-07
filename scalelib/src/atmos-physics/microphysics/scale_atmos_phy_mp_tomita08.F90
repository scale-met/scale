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

  use scale_tracer_tomita08
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_MP_tomita08_setup
  public :: ATMOS_PHY_MP_tomita08
  public :: ATMOS_PHY_MP_tomita08_CloudFraction
  public :: ATMOS_PHY_MP_tomita08_EffectiveRadius
  public :: ATMOS_PHY_MP_tomita08_Mixingratio

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, allocatable :: vterm(:,:,:,:) ! terminal velocity of each tracer [m/s]

  real(RP), public, target :: ATMOS_PHY_MP_DENS(MP_QA) ! hydrometeor density [kg/m3]=[g/L]

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
  logical, private, save  :: MP_doreport_tendency = .false. ! report tendency of each process?
  logical, private, save  :: MP_donegative_fixer  = .true. ! apply negative fixer?

  real(RP), private, save      :: dens00 = 1.28_RP !< standard density [kg/m3]

  ! Parameter for Marshall-Palmer distribution
  real(RP), private, parameter :: N0r = 8.E6_RP !< intercept parameter for rain    [1/m4]
  real(RP), private, parameter :: N0s = 3.E6_RP !< intercept parameter for snow    [1/m4]
  real(RP), private, parameter :: N0g = 4.E6_RP !< intercept parameter for graupel [1/m4]

  real(RP), private, save      :: dens_s = 100.0_RP !< density of snow    [kg/m3]
  real(RP), private, save      :: dens_g = 400.0_RP !< density of graupel [kg/m3]
                                                    !   graupel : 400
                                                    !   hail    : 917
  real(RP), private, save      :: drag_g = 0.6_RP   !< drag coefficient for graupel

  ! Empirical parameter
  real(RP), private, save      :: Ar, As, Ag
  real(RP), private, save      :: Br, Bs, Bg
  real(RP), private, save      :: Cr, Cs, Cg
  real(RP), private, save      :: Dr, Ds, Dg

  ! GAMMA function
  real(RP), private, save      :: GAM, GAM_2, GAM_3

  real(RP), private, save      :: GAM_1br, GAM_2br, GAM_3br
  real(RP), private, save      :: GAM_3dr
  real(RP), private, save      :: GAM_6dr
  real(RP), private, save      :: GAM_1brdr
  real(RP), private, save      :: GAM_5dr_h

  real(RP), private, save      :: GAM_1bs, GAM_2bs, GAM_3bs
  real(RP), private, save      :: GAM_3ds
  real(RP), private, save      :: GAM_1bsds
  real(RP), private, save      :: GAM_5ds_h

  real(RP), private, save      :: GAM_1bg, GAM_3dg
  real(RP), private, save      :: GAM_1bgdg
  real(RP), private, save      :: GAM_5dg_h

  ! Accretion parameter
  real(RP), private, save      :: Eiw = 1.0_RP  !< collection efficiency of cloud ice for cloud water
  real(RP), private, save      :: Erw = 1.0_RP  !< collection efficiency of rain    for cloud water
  real(RP), private, save      :: Esw = 1.0_RP  !< collection efficiency of snow    for cloud water
  real(RP), private, save      :: Egw = 1.0_RP  !< collection efficiency of graupel for cloud water
  real(RP), private, save      :: Eri = 1.0_RP  !< collection efficiency of rain    for cloud ice
  real(RP), private, save      :: Esi = 1.0_RP  !< collection efficiency of snow    for cloud ice
  real(RP), private, save      :: Egi = 0.1_RP  !< collection efficiency of graupel for cloud ice
  real(RP), private, save      :: Esr = 1.0_RP  !< collection efficiency of snow    for rain
  real(RP), private, save      :: Egr = 1.0_RP  !< collection efficiency of graupel for rain
  real(RP), private, save      :: Egs = 1.0_RP  !< collection efficiency of graupel for snow
  real(RP), private, save      :: gamma_sacr = 25.E-3_RP !< effect of low temperature for Esi
  real(RP), private, save      :: gamma_gacs = 90.E-3_RP !< effect of low temperature for Egs

  ! Auto-conversion parameter
  real(RP), private, save      :: Nc_lnd     = 2000.0_RP !< number concentration of cloud water (land)  [1/cc]
  real(RP), private, save      :: Nc_ocn     =   50.0_RP !< number concentration of cloud water (ocean) [1/cc]
  real(RP), private, allocatable :: Nc_def(:,:)          !< number concentration of cloud water         [1/cc]

  real(RP), private, save      :: beta_saut  =  6.E-3_RP  !< auto-conversion factor beta  for ice
  real(RP), private, save      :: gamma_saut = 60.E-3_RP  !< auto-conversion factor gamma for ice
  real(RP), private, save      :: beta_gaut  =  1.E-3_RP  !< auto-conversion factor beta  for snow
  real(RP), private, save      :: gamma_gaut = 90.E-3_RP  !< auto-conversion factor gamma for snow
  real(RP), private, save      :: qicrt_saut =  0.0_RP    !< mixing ratio threshold for Psaut [kg/kg]
  real(RP), private, save      :: qscrt_gaut =  6.E-4_RP  !< mixing ratio threshold for Pgaut [kg/kg]

  ! Evaporation, Sublimation parameter
  real(RP), private, parameter :: Da0    = 2.428E-2_RP !< thermal diffusion coefficient of air at 0C,1atm [J/m/s/K]
  real(RP), private, parameter :: dDa_dT =  7.47E-5_RP !< Coefficient of Da depending on temperature      [J/m/s/K/K]
  real(RP), private, parameter :: Dw0    = 2.222E-5_RP !< diffusion coefficient of water vapor in the air at 0C,1atm [m2/s]
  real(RP), private, parameter :: dDw_dT =  1.37E-7_RP !< Coefficient of Dw depending on temperature                 [m2/s/K]
  real(RP), private, parameter :: mu0    = 1.718E-5_RP !< kinematic viscosity of air at 0C,1atm      [m2/s*kg/m3]
  real(RP), private, parameter :: dmu_dT =  5.28E-8_RP !< Coefficient of mu depending on temperature [m2/s/K*kg/m3]

  real(RP), private, save      :: f1r = 0.78_RP  !< ventilation factor 1 for rain
  real(RP), private, save      :: f2r = 0.27_RP  !< ventilation factor 2 for rain
  real(RP), private, save      :: f1s = 0.65_RP  !< ventilation factor 1 for snow
  real(RP), private, save      :: f2s = 0.39_RP  !< ventilation factor 2 for snow
  real(RP), private, save      :: f1g = 0.78_RP  !< ventilation factor 1 for graupel
  real(RP), private, save      :: f2g = 0.27_RP  !< ventilation factor 2 for graupel

  real(RP), private, save      :: qscrt_sdep =  1.E-12_RP !< mixing ratio threshold for Psdep [kg/kg]
  real(RP), private, save      :: qgcrt_gdep =  1.E-12_RP !< mixing ratio threshold for Pgdep [kg/kg]

  ! Freezing parameter
  real(RP), private, save      :: A_gfrz = 0.66_RP  !< freezing factor [/K]
  real(RP), private, save      :: B_gfrz = 100.0_RP !< freezing factor [/m3/s]

  ! Bergeron process parameter
  real(RP), private, save      :: mi40  = 2.46E-10_RP !< mass              of a 40 micron ice crystal [kg]
  real(RP), private, save      :: mi50  = 4.80E-10_RP !< mass              of a 50 micron ice crystal [kg]
  real(RP), private, save      :: vti50 = 1.0_RP      !< terminal velocity of a 50 micron ice crystal [m/s]
  real(RP), private, save      :: Ri50  = 5.E-5_RP    !< radius            of a 50 micron ice crystal [m]

  integer, private, parameter :: w_nmax = 42
  integer, private, parameter :: I_dqv_dt  =  1 !
  integer, private, parameter :: I_dqc_dt  =  2 !
  integer, private, parameter :: I_dqr_dt  =  3 !
  integer, private, parameter :: I_dqi_dt  =  4 !
  integer, private, parameter :: I_dqs_dt  =  5 !
  integer, private, parameter :: I_dqg_dt  =  6 !
  integer, private, parameter :: I_delta1  =  7 ! separation switch for r->s,g
  integer, private, parameter :: I_delta2  =  8 ! separation switch for s->g
  integer, private, parameter :: I_delta3  =  9 ! separation switch for ice sublimation
  integer, private, parameter :: I_RLMDr   = 10
  integer, private, parameter :: I_RLMDs   = 11
  integer, private, parameter :: I_RLMDg   = 12
  integer, private, parameter :: I_Piacr   = 13 ! r->s,g
  integer, private, parameter :: I_Psacr   = 14 ! r->s,g
  integer, private, parameter :: I_Praci   = 15 ! i->s,g
  integer, private, parameter :: I_Psmlt   = 16 ! s->r
  integer, private, parameter :: I_Pgmlt   = 17 ! g->r
  integer, private, parameter :: I_Praut   = 18 ! c->r
  integer, private, parameter :: I_Pracw   = 19 ! c->r
  integer, private, parameter :: I_Psacw   = 20 ! c->s
  integer, private, parameter :: I_Psfw    = 21 ! c->s
  integer, private, parameter :: I_Pgacw   = 22 ! c->g
  integer, private, parameter :: I_Prevp   = 23 ! r->v
  integer, private, parameter :: I_Piacr_s = 24 ! r->s
  integer, private, parameter :: I_Psacr_s = 25 ! r->s
  integer, private, parameter :: I_Piacr_g = 26 ! r->g
  integer, private, parameter :: I_Psacr_g = 27 ! r->g
  integer, private, parameter :: I_Pgacr   = 28 ! r->g
  integer, private, parameter :: I_Pgfrz   = 29 ! r->g
  integer, private, parameter :: I_Psaut   = 30 ! i->s
  integer, private, parameter :: I_Praci_s = 31 ! i->s
  integer, private, parameter :: I_Psaci   = 32 ! i->s
  integer, private, parameter :: I_Psfi    = 33 ! i->s
  integer, private, parameter :: I_Praci_g = 34 ! i->g
  integer, private, parameter :: I_Pgaci   = 35 ! i->g
  integer, private, parameter :: I_Psdep   = 36 ! v->s
  integer, private, parameter :: I_Pssub   = 37 ! s->v
  integer, private, parameter :: I_Pgaut   = 38 ! s->g
  integer, private, parameter :: I_Pracs   = 39 ! s->g
  integer, private, parameter :: I_Pgacs   = 40 ! s->g
  integer, private, parameter :: I_Pgdep   = 41 ! v->g
  integer, private, parameter :: I_Pgsub   = 42 ! g->v

  real(RP),        private, allocatable :: w(:,:) ! working array
  integer,                private, save :: w_histid (w_nmax) = -999
  logical,                private, save :: w_zinterp(w_nmax) = .false.
  character(len=H_SHORT), private, save :: w_name   (w_nmax)
  character(len=H_MID),   private, save :: w_desc   (w_nmax) = ''
  character(len=H_SHORT), private, save :: w_unit   (w_nmax) = 'kg/kg/s'

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

  real(RP), private, allocatable :: work3D(:,:,:) !< for history output

  logical, private :: debug

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup Cloud Microphysics
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_tomita08_setup( MP_TYPE )
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       PI     => CONST_PI,    &
       GRAV   => CONST_GRAV,  &
       dens_w => CONST_DWATR,  &
       dens_i => CONST_DICE
    use scale_specfunc, only: &
       SF_gamma
    use scale_history, only: &
       HIST_reg
    implicit none
    character(len=H_SHORT), intent(in) :: MP_TYPE

    NAMELIST / PARAM_ATMOS_PHY_MP / &
       MP_doreport_tendency, &
       MP_donegative_fixer

    NAMELIST / PARAM_ATMOS_PHY_MP_TOMITA08 / &
       dens_s, &
       dens_g, &
       debug

    integer :: ierr
    integer :: i, j, ip
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Cloud Microphisics]/Categ[ATMOS]'
    if( IO_L ) write(IO_FID_LOG,*) '*** TOMITA08: 1-moment bulk 6 category'

    if ( MP_TYPE /= 'TOMITA08' ) then
       if ( IO_L ) write(IO_FID_LOG,*) 'xxx ATMOS_PHY_MP_TYPE is not TOMITA08. Check!'
       call PRC_MPIstop
    endif

    allocate( vterm(KA,IA,JA,QA) )
    allocate( Nc_def(IA,JA) )
    allocate( w(w_nmax,KMAX*IMAX*JMAX) )
    allocate( work3D(KA,IA,JA) )


    if (      I_QV <= 0 &
         .OR. I_QC <= 0 &
         .OR. I_QR <= 0 &
         .OR. I_QI <= 0 &
         .OR. I_QS <= 0 &
         .OR. I_QG <= 0 ) then
       if ( IO_L ) write(IO_FID_LOG,*) 'xxx TOMITA08 needs QV/C/R/I/S/G tracer. Check!'
       call PRC_MPIstop
    endif

    !--- empirical coefficients A, B, C, D
    Ar = PI * dens_w / 6.0_RP
    As = PI * dens_s / 6.0_RP
    Ag = PI * dens_g / 6.0_RP

    Br = 3.0_RP
    Bs = 3.0_RP
    Bg = 3.0_RP

    Cr = 130.00_RP
    Cs =   4.84_RP
    Cg = sqrt( ( 4.0_RP * dens_g * GRAV ) / ( 3.0_RP * dens00 * drag_g ) )

    Dr = 0.50_RP
    Ds = 0.25_RP
    Dg = 0.50_RP

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

    vterm(:,:,:,:) = 0.0_RP

    ATMOS_PHY_MP_DENS(I_mp_QC) = dens_w
    ATMOS_PHY_MP_DENS(I_mp_QR) = dens_w
    ATMOS_PHY_MP_DENS(I_mp_QI) = dens_i
    ATMOS_PHY_MP_DENS(I_mp_QS) = dens_s
    ATMOS_PHY_MP_DENS(I_mp_QG) = dens_g

    do j = JS, JE
    do i = IS, IE
       Nc_def(i,j) = Nc_ocn
    enddo
    enddo

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

    do ip = 1, w_nmax
       call HIST_reg( w_histid(ip), w_zinterp(ip), w_name(ip), w_desc(ip), w_unit(ip), ndim=3 )
    enddo

    w_desc(:) = w_name(:)

    work3D(:,:,:) = 0.0_RP

    return
  end subroutine ATMOS_PHY_MP_tomita08_setup

  !-----------------------------------------------------------------------------
  !> Cloud Microphysics
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_tomita08( &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       RHOT, &
       QTRC  )
    use scale_time, only: &
       dt => TIME_DTSEC_ATMOS_PHY_MP
    use scale_history, only: &
       HIST_in
    use scale_atmos_phy_mp_common, only: &
       MP_negative_fixer        => ATMOS_PHY_MP_negative_fixer,       &
       MP_precipitation         => ATMOS_PHY_MP_precipitation,        &
       MP_saturation_adjustment => ATMOS_PHY_MP_saturation_adjustment
    use scale_atmos_thermodyn, only: &
       THERMODYN_rhoe        => ATMOS_THERMODYN_rhoe,       &
       THERMODYN_rhot        => ATMOS_THERMODYN_rhot,       &
       THERMODYN_temp_pres_E => ATMOS_THERMODYN_temp_pres_E
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

    real(RP) :: RHOE_t(KA,IA,JA)
    real(RP) :: QTRC_t(KA,IA,JA,QA)
    real(RP) :: RHOE  (KA,IA,JA)
    real(RP) :: temp  (KA,IA,JA)
    real(RP) :: pres  (KA,IA,JA)

    real(RP) :: flux_tot (KA,IA,JA)
    real(RP) :: flux_rain(KA,IA,JA)
    real(RP) :: flux_snow(KA,IA,JA)
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Microphysics(tomita08)'

    if ( MP_donegative_fixer ) then
       call MP_negative_fixer( DENS(:,:,:),  & ! [INOUT]
                               RHOT(:,:,:),  & ! [INOUT]
                               QTRC(:,:,:,:) ) ! [INOUT]
    endif

    call THERMODYN_rhoe( RHOE(:,:,:),  & ! [OUT]
                         RHOT(:,:,:),  & ! [IN]
                         QTRC(:,:,:,:) ) ! [IN]

    !##### MP Main #####
    RHOE_t(:,:,:)   = 0.0_RP
    QTRC_t(:,:,:,:) = 0.0_RP

    call MP_tomita08( RHOE_t(:,:,:),   & ! [OUT]
                      QTRC_t(:,:,:,:), & ! [OUT]
                      RHOE  (:,:,:),   & ! [INOUT]
                      QTRC  (:,:,:,:), & ! [INOUT]
                      DENS  (:,:,:)    ) ! [IN]

    call MP_tomita08_vterm( vterm(:,:,:,:), & ! [OUT]
                            DENS (:,:,:),   & ! [IN]
                            QTRC (:,:,:,:)  ) ! [IN]

    call THERMODYN_temp_pres_E( temp(:,:,:),  & ! [OUT]
                                pres(:,:,:),  & ! [OUT]
                                DENS(:,:,:),  & ! [IN]
                                RHOE(:,:,:),  & ! [IN]
                                QTRC(:,:,:,:) ) ! [IN]

    call MP_precipitation( flux_rain(:,:,:), & ! [OUT]
                           flux_snow(:,:,:), & ! [OUT]
                           DENS (:,:,:),     & ! [INOUT]
                           MOMZ (:,:,:),     & ! [INOUT]
                           MOMX (:,:,:),     & ! [INOUT]
                           MOMY (:,:,:),     & ! [INOUT]
                           RHOE (:,:,:),     & ! [INOUT]
                           QTRC (:,:,:,:),   & ! [INOUT]
                           vterm(:,:,:,:),   & ! [IN]
                           temp (:,:,:),     & ! [IN]
                           dt                ) ! [IN]

    call MP_saturation_adjustment( RHOE_t(:,:,:),   & ! [INOUT]
                                   QTRC_t(:,:,:,:), & ! [INOUT]
                                   RHOE(:,:,:),     & ! [INOUT]
                                   QTRC(:,:,:,:),   & ! [INOUT]
                                   DENS(:,:,:)      ) ! [IN]

    if ( MP_doreport_tendency ) then
       call HIST_in( QTRC_t(:,:,:,I_QV), 'QV_t_mp', 'tendency of QV in mp', 'kg/kg/s', dt )
       call HIST_in( QTRC_t(:,:,:,I_QC), 'QC_t_mp', 'tendency of QC in mp', 'kg/kg/s', dt )
       call HIST_in( QTRC_t(:,:,:,I_QR), 'QR_t_mp', 'tendency of QR in mp', 'kg/kg/s', dt )
       call HIST_in( QTRC_t(:,:,:,I_QI), 'QI_t_mp', 'tendency of QI in mp', 'kg/kg/s', dt )
       call HIST_in( QTRC_t(:,:,:,I_QS), 'QS_t_mp', 'tendency of QS in mp', 'kg/kg/s', dt )
       call HIST_in( QTRC_t(:,:,:,I_QG), 'QG_t_mp', 'tendency of QG in mp', 'kg/kg/s', dt )

       call HIST_in( RHOE_t(:,:,:), 'RHOE_t_mp', 'tendency of rhoe in mp', 'J/m3/s', dt )

       call HIST_in( vterm(:,:,:,I_QR), 'Vterm_QR', 'terminal velocity of QR', 'm/s', dt )
       call HIST_in( vterm(:,:,:,I_QI), 'Vterm_QI', 'terminal velocity of QI', 'm/s', dt )
       call HIST_in( vterm(:,:,:,I_QS), 'Vterm_QS', 'terminal velocity of QS', 'm/s', dt )
       call HIST_in( vterm(:,:,:,I_QG), 'Vterm_QG', 'terminal velocity of QG', 'm/s', dt )
    endif

    !##### END MP Main #####

    call THERMODYN_rhot( RHOT(:,:,:),  & ! [OUT]
                         RHOE(:,:,:),  & ! [IN]
                         QTRC(:,:,:,:) ) ! [IN]

    if ( MP_donegative_fixer ) then
       call MP_negative_fixer( DENS(:,:,:),  & ! [INOUT]
                               RHOT(:,:,:),  & ! [INOUT]
                               QTRC(:,:,:,:) ) ! [INOUT]
    endif

    flux_tot(:,:,:) = flux_rain(:,:,:) + flux_snow(:,:,:)
    call HIST_in( flux_rain(KS-1,:,:), 'RAIN', 'surface rain rate', 'kg/m2/s', dt)
    call HIST_in( flux_snow(KS-1,:,:), 'SNOW', 'surface snow rate', 'kg/m2/s', dt)
    call HIST_in( flux_tot (KS-1,:,:), 'PREC', 'surface precipitation rate', 'kg/m2/s', dt)

    return
  end subroutine ATMOS_PHY_MP_tomita08

  !-----------------------------------------------------------------------------
  !> Lin-type cold rain microphysics
  !-----------------------------------------------------------------------------
  subroutine MP_tomita08( &
       RHOE_t, &
       QTRC_t, &
       RHOE0,  &
       QTRC0,  &
       DENS0   )
    use scale_const, only: &
       PI    => CONST_PI,    &
       EPS   => CONST_EPS,   &
       Rvap  => CONST_Rvap,  &
       CL    => CONST_CL,    &
       LHV00 => CONST_LH00,  &
       LHF00 => CONST_LHF00, &
       LHV0  => CONST_LH0,   &
       LHS0  => CONST_LHS0,  &
       LHF0  => CONST_LHF0,  &
       TEM00 => CONST_TEM00, &
       PRE00 => CONST_PRE00
    use scale_time, only: &
       dt => TIME_DTSEC_ATMOS_PHY_MP
    use scale_history, only: &
       HIST_put, &
       HIST_in
    use scale_atmos_thermodyn, only: &
       THERMODYN_temp_pres_E => ATMOS_THERMODYN_temp_pres_E
    use scale_atmos_saturation, only: &
       SATURATION_dens2qsat_liq => ATMOS_SATURATION_dens2qsat_liq, &
       SATURATION_dens2qsat_ice => ATMOS_SATURATION_dens2qsat_ice
    implicit none

    real(RP), intent(inout) :: RHOE_t(KA,IA,JA)    ! tendency rhoe             [J/m3/s]
    real(RP), intent(inout) :: QTRC_t(KA,IA,JA,QA) ! tendency tracer           [kg/kg/s]
    real(RP), intent(inout) :: RHOE0 (KA,IA,JA)    ! density * internal energy [J/m3]
    real(RP), intent(inout) :: QTRC0 (KA,IA,JA,QA) ! mass concentration        [kg/kg]
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
    real(RP) :: q(QA)
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

    real(RP) :: Vti, Vtr, Vts, Vtg  !< terminal velocity
    real(RP) :: Esi_mod, Egs_mod    !< modified accretion efficiency
    real(RP) :: rhoqc               !< rho * qc
    real(RP) :: Nc(KA,IA,JA)        !< Number concentration of cloud water [1/cc]
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

    integer :: ijk_warm, ijk_cold
    integer :: index_warm(KMAX*IMAX*JMAX)
    integer :: index_cold(KMAX*IMAX*JMAX)

    real(RP) :: net, fac, fac_sw
    real(RP) :: tend(I_QV:I_QG,KMAX*IMAX*JMAX)

    real(RP) :: zerosw
    real(RP) :: tmp

    integer :: k, i, j, iq, ijk, indirect, ip
    !---------------------------------------------------------------------------

    call PROF_rapstart('MP_tomita08')

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       Nc(k,i,j) = Nc_def(i,j)
    enddo
    enddo
    enddo

    rdt = 1.0_RP / dt

    call THERMODYN_temp_pres_E( TEMP0(:,:,:),  & ! [OUT]
                                PRES0(:,:,:),  & ! [OUT]
                                DENS0(:,:,:),  & ! [IN]
                                RHOE0(:,:,:),  & ! [IN]
                                QTRC0(:,:,:,:) ) ! [IN]

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

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       ijk = ( j - JS ) * KMAX * IMAX &
           + ( i - IS ) * KMAX        &
           + ( k - KS )               &
           + 1

       ! store to work
       dens = DENS0(k,i,j)
       rhoe = RHOE0(k,i,j)
       temp = TEMP0(k,i,j)
       pres = PRES0(k,i,j)
       do iq = I_QV, I_QG
          q(iq) = QTRC0(k,i,j,iq)
       enddo

       ! saturation ratio S
       Sliq = q(I_QV) / max( QSATL(k,i,j), EPS )
       Sice = q(I_QV) / max( QSATI(k,i,j), EPS )

       Rdens    = 1.0_RP / dens
       rho_fact = sqrt( dens00 * Rdens )
       temc     = temp - TEM00

       w(I_delta1,ijk) = ( 0.5_RP + sign(0.5_RP, q(I_QR) - 1.E-4_RP ) )

       w(I_delta2,ijk) = ( 0.5_RP + sign(0.5_RP, 1.E-4_RP - q(I_QR) ) ) &
                       * ( 0.5_RP + sign(0.5_RP, 1.E-4_RP - q(I_QS) ) )

       w(I_delta3,ijk) = 0.5_RP + sign(0.5_RP, Sice - 1.0_RP )

       w(I_dqv_dt,ijk) = q(I_QV) * rdt
       w(I_dqc_dt,ijk) = q(I_QC) * rdt
       w(I_dqr_dt,ijk) = q(I_QR) * rdt
       w(I_dqi_dt,ijk) = q(I_QI) * rdt
       w(I_dqs_dt,ijk) = q(I_QS) * rdt
       w(I_dqg_dt,ijk) = q(I_QG) * rdt

       Bergeron_sw = ( 0.5_RP + sign(0.5_RP, temc + 30.0_RP ) ) &
                   * ( 0.5_RP + sign(0.5_RP, 0.0_RP - temc  ) )

       ! slope parameter lambda
       zerosw = 0.5_RP - sign(0.5_RP, q(I_QR) - 1.E-12_RP )
       RLMDr = sqrt(sqrt( dens * q(I_QR) / ( Ar * N0r * GAM_1br ) + zerosw )) * ( 1.0_RP - zerosw )

       zerosw = 0.5_RP - sign(0.5_RP, q(I_QS) - 1.E-12_RP )
       RLMDs = sqrt(sqrt( dens * q(I_QS) / ( As * N0s * GAM_1bs ) + zerosw )) * ( 1.0_RP - zerosw )

       zerosw = 0.5_RP - sign(0.5_RP, q(I_QG) - 1.E-12_RP )
       RLMDg = sqrt(sqrt( dens * q(I_QG) / ( Ag * N0g * GAM_1bg ) + zerosw )) * ( 1.0_RP - zerosw )

       w(I_RLMDr,ijk) = RLMDr
       w(I_RLMDs,ijk) = RLMDs
       w(I_RLMDg,ijk) = RLMDg

       RLMDr_dr = sqrt( RLMDr )       ! **Dr
       RLMDs_ds = sqrt( sqrt(RLMDs) ) ! **Ds
       RLMDg_dg = sqrt( RLMDg )       ! **Dg

       RLMDr_2   = RLMDr**2
       RLMDr_3   = RLMDr**3
       RLMDr_7   = RLMDr**7
       RLMDr_1br = RLMDr**4 ! (1+Br)
       RLMDr_2br = RLMDr**5 ! (2+Br)
       RLMDr_3br = RLMDr**6 ! (3+Br)
       RLMDr_3dr = RLMDr**3 * RLMDr_dr
       RLMDr_5dr = RLMDr**5 * RLMDr_dr
       RLMDr_6dr = RLMDr**6 * RLMDr_dr

       RLMDs_2   = RLMDs**2
       RLMDs_3   = RLMDs**3
       RLMDs_1bs = RLMDs**4 ! (1+Bs)
       RLMDs_2bs = RLMDs**5 ! (2+Bs)
       RLMDs_3bs = RLMDs**6 ! (3+Bs)
       RLMDs_3ds = RLMDs**3 * RLMDs_ds
       RLMDs_5ds = RLMDs**5 * RLMDs_ds

       RLMDg_2   = RLMDg**2
       RLMDg_3   = RLMDg**3
       RLMDg_3dg = RLMDg**3 * RLMDg_dg
       RLMDg_5dg = RLMDg**5 * RLMDg_dg

       !---< terminal velocity >
       zerosw = 0.5_RP + sign(0.5_RP, q(I_QI) - 1.E-8_RP )
       Vti = -3.29_RP * ( dens * q(I_QI) * zerosw )**0.16_RP
       Vtr = -Cr * rho_fact * GAM_1brdr / GAM_1br * RLMDr_dr
       Vts = -Cs * rho_fact * GAM_1bsds / GAM_1bs * RLMDs_ds
       Vtg = -Cg * rho_fact * GAM_1bgdg / GAM_1bg * RLMDg_dg

       !---< Accretion >---
       Esi_mod = min( Esi, Esi * exp( gamma_sacr * temc ) )
       Egs_mod = min( Egs, Egs * exp( gamma_gacs * temc ) )

       ! [Pracw] accretion rate of cloud water by rain
       w(I_Pracw,ijk) = q(I_QC) * 0.25_RP * PI * Erw * N0r * Cr * GAM_3dr * RLMDr_3dr * rho_fact

       ! [Psacw] accretion rate of cloud water by snow
       w(I_Psacw,ijk) = q(I_QC) * 0.25_RP * PI * Esw * N0s * Cs * GAM_3ds * RLMDs_3ds * rho_fact

       ! [Pgacw] accretion rate of cloud water by graupel
       w(I_Pgacw,ijk) = q(I_QC) * 0.25_RP * PI * Egw * N0g * Cg * GAM_3dg * RLMDg_3dg * rho_fact

       ! [Praci] accretion rate of cloud ice by rain
       w(I_Praci,ijk) = q(I_QI) * 0.25_RP * PI * Eri * N0r * Cr * GAM_3dr * RLMDr_3dr * rho_fact

       ! [Psaci] accretion rate of cloud ice by snow
       w(I_Psaci,ijk) = q(I_QI) * 0.25_RP * PI * Esi_mod * N0s * Cs * GAM_3ds * RLMDs_3ds * rho_fact

       ! [Pgaci] accretion rate of cloud ice by grupel
       w(I_Pgaci,ijk) = q(I_QI) * 0.25_RP * PI * Egi * N0g * Cg * GAM_3dg * RLMDg_3dg * rho_fact

       ! [Piacr] accretion rate of rain by cloud ice
       w(I_Piacr,ijk) = q(I_QI) * Ar / 4.19E-13_RP * 0.25_RP * PI * Eri * N0r * Cr * GAM_6dr * RLMDr_6dr * rho_fact

       ! [Psacr] accretion rate of rain by snow
       w(I_Psacr,ijk) = Ar * 0.25_RP * PI * Rdens * Esr * N0r * N0s * abs(Vtr-Vts) &
                      * (          GAM_1br * GAM_3 * RLMDr_1br * RLMDs_3 &
                        + 2.0_RP * GAM_2br * GAM_2 * RLMDr_2br * RLMDs_2 &
                        +          GAM_3br * GAM   * RLMDr_3br * RLMDs   )

       ! [Pgacr] accretion rate of rain by graupel
       w(I_Pgacr,ijk) = Ar * 0.25_RP * PI * Rdens * Egr * N0g * N0r * abs(Vtg-Vtr) &
                      * (          GAM_1br * GAM_3 * RLMDr_1br * RLMDg_3 &
                        + 2.0_RP * GAM_2br * GAM_2 * RLMDr_2br * RLMDg_2 &
                        +          GAM_3br * GAM   * RLMDr_3br * RLMDg   )

       ! [Pracs] accretion rate of snow by rain
       w(I_Pracs,ijk) = As * 0.25_RP * PI * Rdens * Esr * N0s * N0r * abs(Vtr-Vts) &
                      * (          GAM_1bs * GAM_3 * RLMDs_1bs * RLMDr_3 &
                        + 2.0_RP * GAM_2bs * GAM_2 * RLMDs_2bs * RLMDr_2 &
                        +          GAM_3bs * GAM   * RLMDs_3bs * RLMDr   )

       ! [Pgacs] accretion rate of snow by graupel
       w(I_Pgacs,ijk) = As * 0.25_RP * PI * Rdens * Egs_mod * N0g * N0s * abs(Vtg-Vts) &
                      * (          GAM_1bs * GAM_3 * RLMDs_1bs * RLMDg_3 &
                        + 2.0_RP * GAM_2bs * GAM_2 * RLMDs_2bs * RLMDg_2 &
                        +          GAM_3bs * GAM   * RLMDs_3bs * RLMDg   )

       !---< Auto-conversion >---
       ! [Praut] auto-conversion rate from cloud water to rain
       rhoqc = dens * q(I_QC) * 1000.0_RP ! [g/m3]
       Dc    = 0.146_RP - 5.964E-2_RP * log( Nc(k,i,j) / 2000.0_RP )
       w(I_Praut,ijk) = Rdens * 1.67E-5_RP * rhoqc * rhoqc / ( 5.0_RP + 3.66E-2_RP * Nc(k,i,j) / ( Dc * rhoqc + EPS ) )

       ! [Psaut] auto-conversion rate from cloud ice to snow
       betai = min( beta_saut, beta_saut * exp( gamma_saut * temc ) )
       w(I_Psaut,ijk) = max( betai*(q(I_QI)-qicrt_saut), 0.0_RP )

       ! [Pgaut] auto-conversion rate from snow to graupel
       betas = min( beta_gaut, beta_gaut * exp( gamma_gaut * temc ) )
       w(I_Pgaut,ijk) = max( betas*(q(I_QS)-qscrt_gaut), 0.0_RP )

       !---< Evaporation, Sublimation >---
       Da = ( Da0 + dDa_dT * temc )
       Kd = ( Dw0 + dDw_dT * temc ) * PRE00 / pres
       NU = ( mu0 + dmu_dT * temc ) * Rdens

       Glv  = 1.0_RP / ( LHV0/(Da*temp) * ( LHV0/(Rvap*temp) - 1.0_RP ) + 1.0_RP/(Kd*dens*QSATL(k,i,j)) )
       Giv  = 1.0_RP / ( LHS0/(Da*temp) * ( LHS0/(Rvap*temp) - 1.0_RP ) + 1.0_RP/(Kd*dens*QSATI(k,i,j)) )
       Gil  = 1.0_RP / LHF0 * (Da*temc)

       ! [Prevp] evaporation rate of rain
       ventr = f1r * GAM_2 * RLMDr_2 + f2r * sqrt( Cr * rho_fact / NU * RLMDr_5dr ) * GAM_5dr_h

       w(I_Prevp,ijk) = 2.0_RP * PI * Rdens * N0r * ( 1.0_RP-min(Sliq,1.0_RP) ) * Glv * ventr

       ! [Psdep,Pssub] deposition/sublimation rate for snow
       vents = f1s * GAM_2 * RLMDs_2 + f2s * sqrt( Cs * rho_fact / NU * RLMDs_5ds ) * GAM_5ds_h

       tmp = 2.0_RP * PI * Rdens * N0s * ( Sice-1.0_RP ) * Giv * vents

       w(I_Psdep,ijk) = ( w(I_delta3,ijk)        ) * tmp ! Sice < 1
       w(I_Pssub,ijk) = ( w(I_delta3,ijk)-1.0_RP ) * tmp ! Sice > 1

       ! [Psmlt] melting rate of snow
       w(I_Psmlt,ijk) = 2.0_RP * PI * Rdens * N0s * Gil * vents &
                      + CL * temc / LHF0 * ( w(I_Psacw,ijk) + w(I_Psacr,ijk) )

       ! [Pgdep/pgsub] deposition/sublimation rate for graupel
       ventg = f1g * GAM_2 * RLMDg_2 + f2g * sqrt( Cg * rho_fact / NU * RLMDg_5dg ) * GAM_5dg_h

       tmp = 2.0_RP * PI * Rdens * N0g * ( Sice-1.0_RP ) * Giv * ventg

       w(I_Pgdep,ijk) = ( w(I_delta3,ijk)        ) * tmp ! Sice < 1
       w(I_Pgsub,ijk) = ( w(I_delta3,ijk)-1.0_RP ) * tmp ! Sice > 1

       ! [Pgmlt] melting rate of graupel
       w(I_Pgmlt,ijk) = 2.0_RP * PI * Rdens * N0g * Gil * ventg &
                      + CL * temc / LHF0 * ( w(I_Pgacw,ijk) + w(I_Pgacr,ijk) )

       ! [Pgfrz] freezing rate of graupel
       w(I_Pgfrz,ijk) = 2.0_RP * PI * Rdens * N0r * 60.0_RP * B_gfrz * Ar * ( exp(-A_gfrz*temc) - 1.0_RP ) * RLMDr_7

       ! [Psfw,Psfi] ( Bergeron process ) growth rate of snow by Bergeron process from cloud water/ice
       dt1  = ( mi50**ma2(k,i,j) - mi40**ma2(k,i,j) ) / ( a1(k,i,j) * ma2(k,i,j) )
       Ni50 = q(I_QI) * dt / ( mi50 * dt1 )

       w(I_Psfw,ijk) = Bergeron_sw * Ni50 * ( a1(k,i,j) * mi50**a2(k,i,j) &
                                            + PI * Eiw * dens * q(I_QC) * Ri50*Ri50 * vti50 )
       w(I_Psfi,ijk) = Bergeron_sw * q(I_QI) / dt1

    enddo
    enddo
    enddo

!    do ip = 1, w_nmax
!       if( IO_L ) write(IO_FID_LOG,*) w_name(ip), "MAX/MIN:", maxval(w(ip,:)), minval(w(ip,:))
!    enddo

    do ijk = 1, KMAX*IMAX*JMAX
       w(I_Psdep,ijk) = min( w(I_Psdep,ijk), w(I_dqv_dt,ijk) )
       w(I_Pgdep,ijk) = min( w(I_Pgdep,ijk), w(I_dqv_dt,ijk) )

       w(I_Praut,ijk) = min( w(I_Praut,ijk), w(I_dqc_dt,ijk) )
       w(I_Pracw,ijk) = min( w(I_Pracw,ijk), w(I_dqc_dt,ijk) )
       w(I_Psacw,ijk) = min( w(I_Psacw,ijk), w(I_dqc_dt,ijk) )
       w(I_Pgacw,ijk) = min( w(I_Pgacw,ijk), w(I_dqc_dt,ijk) )
       w(I_Psfw ,ijk) = min( w(I_Psfw ,ijk), w(I_dqc_dt,ijk) )

       w(I_Prevp,ijk) = min( w(I_Prevp,ijk), w(I_dqr_dt,ijk) )
       w(I_Piacr,ijk) = min( w(I_Piacr,ijk), w(I_dqr_dt,ijk) )
       w(I_Psacr,ijk) = min( w(I_Psacr,ijk), w(I_dqr_dt,ijk) )
       w(I_Pgacr,ijk) = min( w(I_Pgacr,ijk), w(I_dqr_dt,ijk) )
       w(I_Pgfrz,ijk) = min( w(I_Pgfrz,ijk), w(I_dqr_dt,ijk) )

       w(I_Psaut,ijk) = min( w(I_Psaut,ijk), w(I_dqi_dt,ijk) )
       w(I_Praci,ijk) = min( w(I_Praci,ijk), w(I_dqi_dt,ijk) )
       w(I_Psaci,ijk) = min( w(I_Psaci,ijk), w(I_dqi_dt,ijk) )
       w(I_Pgaci,ijk) = min( w(I_Pgaci,ijk), w(I_dqi_dt,ijk) )
       w(I_Psfi ,ijk) = min( w(I_Psfi ,ijk), w(I_dqi_dt,ijk) )

       w(I_Pgaut,ijk) = min( w(I_Pgaut,ijk), w(I_dqs_dt,ijk) )
       w(I_Pracs,ijk) = min( w(I_Pracs,ijk), w(I_dqs_dt,ijk) )
       w(I_Pgacs,ijk) = min( w(I_Pgacs,ijk), w(I_dqs_dt,ijk) )
       w(I_Psmlt,ijk) = max( w(I_Psmlt,ijk), 0.0_RP          )
       w(I_Psmlt,ijk) = min( w(I_Psmlt,ijk), w(I_dqs_dt,ijk) )
       w(I_Pssub,ijk) = min( w(I_Pssub,ijk), w(I_dqs_dt,ijk) )

       w(I_Pgmlt,ijk) = max( w(I_Pgmlt,ijk), 0.0_RP          )
       w(I_Pgmlt,ijk) = min( w(I_Pgmlt,ijk), w(I_dqg_dt,ijk) )
       w(I_Pgsub,ijk) = min( w(I_Pgsub,ijk), w(I_dqg_dt,ijk) )

       w(I_Piacr_s,ijk) = ( 1.0_RP - w(I_delta1,ijk) ) * w(I_Piacr,ijk)
       w(I_Piacr_g,ijk) = (          w(I_delta1,ijk) ) * w(I_Piacr,ijk)
       w(I_Praci_s,ijk) = ( 1.0_RP - w(I_delta1,ijk) ) * w(I_Praci,ijk)
       w(I_Praci_g,ijk) = (          w(I_delta1,ijk) ) * w(I_Praci,ijk)
       w(I_Psacr_s,ijk) = (          w(I_delta2,ijk) ) * w(I_Psacr,ijk)
       w(I_Psacr_g,ijk) = ( 1.0_RP - w(I_delta2,ijk) ) * w(I_Psacr,ijk)
       w(I_Pracs  ,ijk) = ( 1.0_RP - w(I_delta2,ijk) ) * w(I_Pracs,ijk)
    enddo

    ijk_warm = 0
    ijk_cold = 0
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       ijk = ( j - JS ) * KMAX * IMAX &
           + ( i - IS ) * KMAX        &
           + ( k - KS )               &
           + 1

       if ( TEMP0(k,i,j)-TEM00 > 0.0_RP ) then ! warm rain
          ijk_warm = ijk_warm + 1
          index_warm(ijk_warm) = ijk
       else
          ijk_cold = ijk_cold + 1
          index_cold(ijk_cold) = ijk
       endif
    enddo
    enddo
    enddo

!    if( IO_L ) write(IO_FID_LOG,*) "ijk_warm/cold:", ijk_warm, ijk_cold

    !---< Solve tendencies (Warm Rain) >---
    do indirect = 1, ijk_warm
       ijk = index_warm(indirect)

       ! [QC]
       net = &                ! [prod]
           - w(I_Praut,ijk) & ! [loss] c->r
           - w(I_Pracw,ijk) & ! [loss] c->r
           - w(I_Psacw,ijk) & ! [loss] c->s=r
           - w(I_Pgacw,ijk)   ! [loss] c->g=r

       fac_sw = 0.5_RP + sign( 0.5_RP, net+EPS ) ! if production > loss , fac_sw=1
       fac    = (          fac_sw ) &
              + ( 1.0_RP - fac_sw ) * min( -w(I_dqc_dt,ijk)/(net-fac_sw), 1.0_RP ) ! loss limiter

       w(I_Praut,ijk) = w(I_Praut,ijk) * fac
       w(I_Pracw,ijk) = w(I_Pracw,ijk) * fac
       w(I_Psacw,ijk) = w(I_Psacw,ijk) * fac
       w(I_Pgacw,ijk) = w(I_Pgacw,ijk) * fac

       ! [QI]
       ! no [prod] & [loss]

       ! [QR]
       net = w(I_Praut,ijk) & ! [prod] c->r
           + w(I_Pracw,ijk) & ! [prod] c->r
           + w(I_Psacw,ijk) & ! [prod] c->s=r
           + w(I_Pgacw,ijk) & ! [prod] c->g=r
           + w(I_Psmlt,ijk) & ! [prod] s->r
           + w(I_Pgmlt,ijk) & ! [prod] g->r
           - w(I_Prevp,ijk)   ! [loss] r->v

       fac_sw = 0.5_RP + sign( 0.5_RP, net+EPS ) ! if production > loss , fac_sw=1
       fac    = (          fac_sw ) &
              + ( 1.0_RP - fac_sw ) * min( -w(I_dqr_dt,ijk)/(net-fac_sw), 1.0_RP ) ! loss limiter

       w(I_Praut,ijk) = w(I_Praut,ijk) * fac
       w(I_Pracw,ijk) = w(I_Pracw,ijk) * fac
       w(I_Psacw,ijk) = w(I_Psacw,ijk) * fac
       w(I_Pgacw,ijk) = w(I_Pgacw,ijk) * fac
       w(I_Psmlt,ijk) = w(I_Psmlt,ijk) * fac
       w(I_Pgmlt,ijk) = w(I_Pgmlt,ijk) * fac
       w(I_Prevp,ijk) = w(I_Prevp,ijk) * fac

       ! [QS]
       net = &                ! [prod]
           - w(I_Psmlt,ijk) & ! [loss] s->r
           - w(I_Pgacs,ijk)   ! [loss] s->g

       fac_sw = 0.5_RP + sign( 0.5_RP, net+EPS ) ! if production > loss , fac_sw=1
       fac    = (          fac_sw ) &
              + ( 1.0_RP - fac_sw ) * min( -w(I_dqs_dt,ijk)/(net-fac_sw), 1.0_RP ) ! loss limiter

       w(I_Psmlt,ijk) = w(I_Psmlt,ijk) * fac
       w(I_Pgacs,ijk) = w(I_Pgacs,ijk) * fac

       ! [QG]
       net = w(I_Pgacs,ijk) & ! [prod] s->g
           - w(I_Pgmlt,ijk)   ! [loss] g->r

       fac_sw = 0.5_RP + sign( 0.5_RP, net+EPS ) ! if production > loss , fac_sw=1
       fac    = (          fac_sw ) &
              + ( 1.0_RP - fac_sw ) * min( -w(I_dqg_dt,ijk)/(net-fac_sw), 1.0_RP ) ! loss limiter

       w(I_Pgacs,ijk) = w(I_Pgacs,ijk) * fac
       w(I_Pgmlt,ijk) = w(I_Pgmlt,ijk) * fac

       net = &                ! [prod]
           - w(I_Psmlt,ijk) & ! [loss] s->r
           - w(I_Pgacs,ijk)   ! [loss] s->g
    enddo



    do indirect = 1, ijk_warm
       ijk = index_warm(indirect)

       ! re-calc net production & loss
       tend(I_QC,ijk) = &                ! [prod]
                      - w(I_Praut,ijk) & ! [loss] c->r
                      - w(I_Pracw,ijk) & ! [loss] c->r
                      - w(I_Psacw,ijk) & ! [loss] c->s=r
                      - w(I_Pgacw,ijk)   ! [loss] c->g=r
       tend(I_QR,ijk) = w(I_Praut,ijk) & ! [prod] c->r
                      + w(I_Pracw,ijk) & ! [prod] c->r
                      + w(I_Psacw,ijk) & ! [prod] c->s=r
                      + w(I_Pgacw,ijk) & ! [prod] c->g=r
                      + w(I_Psmlt,ijk) & ! [prod] s->r
                      + w(I_Pgmlt,ijk) & ! [prod] g->r
                      - w(I_Prevp,ijk)   ! [loss] r->v
       tend(I_QI,ijk) = 0.0_RP
       tend(I_QS,ijk) = &                ! [prod]
                      - w(I_Psmlt,ijk) & ! [loss] s->r
                      - w(I_Pgacs,ijk)   ! [loss] s->g
       tend(I_QG,ijk) = w(I_Pgacs,ijk) & ! [prod] s->g
                      - w(I_Pgmlt,ijk)   ! [loss] g->r

       tend(I_QC,ijk) = max( tend(I_QC,ijk), -w(I_dqc_dt,ijk) )
       tend(I_QR,ijk) = max( tend(I_QR,ijk), -w(I_dqr_dt,ijk) )
       tend(I_QI,ijk) = max( tend(I_QI,ijk), -w(I_dqi_dt,ijk) )
       tend(I_QS,ijk) = max( tend(I_QS,ijk), -w(I_dqs_dt,ijk) )
       tend(I_QG,ijk) = max( tend(I_QG,ijk), -w(I_dqg_dt,ijk) )

       tend(I_QV,ijk) = - ( tend(I_QC,ijk) &
                          + tend(I_QR,ijk) &
                          + tend(I_QI,ijk) &
                          + tend(I_QS,ijk) &
                          + tend(I_QG,ijk) )
    enddo

    !---< Solve tendencies (Cold Rain) >---
    do indirect = 1, ijk_cold
       ijk = index_cold(indirect)

       ! [QC]
       net = &                  ! [prod]
           - w(I_Praut  ,ijk) & ! [loss] c->r
           - w(I_Pracw  ,ijk) & ! [loss] c->r
           - w(I_Psacw  ,ijk) & ! [loss] c->s
           - w(I_Psfw   ,ijk) & ! [loss] c->s
           - w(I_Pgacw  ,ijk)   ! [loss] c->g

       fac_sw = 0.5_RP + sign( 0.5_RP, net+EPS ) ! if production > loss , fac_sw=1
       fac    = (          fac_sw ) &
              + ( 1.0_RP - fac_sw ) * min( -w(I_dqc_dt,ijk)/(net-fac_sw), 1.0_RP ) ! loss limiter

       w(I_Praut,ijk) = w(I_Praut,ijk) * fac
       w(I_Pracw,ijk) = w(I_Pracw,ijk) * fac
       w(I_Psacw,ijk) = w(I_Psacw,ijk) * fac
       w(I_Psfw ,ijk) = w(I_Psfw ,ijk) * fac
       w(I_Pgacw,ijk) = w(I_Pgacw,ijk) * fac

       ! [QI]
       net = &                  ! [prod]
           - w(I_Psaut  ,ijk) & ! [loss] i->s
           - w(I_Praci_s,ijk) & ! [loss] i->s
           - w(I_Psaci  ,ijk) & ! [loss] i->s
           - w(I_Psfi   ,ijk) & ! [loss] i->s
           - w(I_Praci_g,ijk) & ! [loss] i->g
           - w(I_Pgaci  ,ijk)   ! [loss] i->g

       fac_sw = 0.5_RP + sign( 0.5_RP, net+EPS ) ! if production > loss , fac_sw=1
       fac    = (          fac_sw ) &
              + ( 1.0_RP - fac_sw ) * min( -w(I_dqi_dt,ijk)/(net-fac_sw), 1.0_RP ) ! loss limiter

       w(I_Praci_s,ijk) = w(I_Praci_s,ijk) * fac
       w(I_Psaut  ,ijk) = w(I_Psaut  ,ijk) * fac
       w(I_Psaci  ,ijk) = w(I_Psaci  ,ijk) * fac
       w(I_Psfi   ,ijk) = w(I_Psfi   ,ijk) * fac
       w(I_Praci_g,ijk) = w(I_Praci_g,ijk) * fac
       w(I_Pgaci  ,ijk) = w(I_Pgaci  ,ijk) * fac

       ! [QR]
       net = w(I_Praut  ,ijk) & ! [prod] c->r
           + w(I_Pracw  ,ijk) & ! [prod] c->r
           - w(I_Prevp  ,ijk) & ! [loss] r->v
           - w(I_Piacr_s,ijk) & ! [loss] r->s
           - w(I_Psacr_s,ijk) & ! [loss] r->s
           - w(I_Piacr_g,ijk) & ! [loss] r->g
           - w(I_Psacr_g,ijk) & ! [loss] r->g
           - w(I_Pgacr  ,ijk) & ! [loss] r->g
           - w(I_Pgfrz  ,ijk)   ! [loss] r->g

       fac_sw = 0.5_RP + sign( 0.5_RP, net+EPS ) ! if production > loss , fac_sw=1
       fac    = (          fac_sw ) &
              + ( 1.0_RP - fac_sw ) * min( -w(I_dqr_dt,ijk)/(net-fac_sw), 1.0_RP ) ! loss limiter

       w(I_Praut  ,ijk) = w(I_Praut  ,ijk) * fac
       w(I_Pracw  ,ijk) = w(I_Pracw  ,ijk) * fac
       w(I_Prevp  ,ijk) = w(I_Prevp  ,ijk) * fac
       w(I_Piacr_s,ijk) = w(I_Piacr_s,ijk) * fac
       w(I_Psacr_s,ijk) = w(I_Psacr_s,ijk) * fac
       w(I_Piacr_g,ijk) = w(I_Piacr_g,ijk) * fac
       w(I_Psacr_g,ijk) = w(I_Psacr_g,ijk) * fac
       w(I_Pgacr  ,ijk) = w(I_Pgacr  ,ijk) * fac
       w(I_Pgfrz  ,ijk) = w(I_Pgfrz  ,ijk) * fac

       !### in the case of Sice-1 > 0, limiter is applied to qv before qs,qg
       ! [QV]
       net = w(I_Prevp,ijk) & ! [prod] r->v
           + w(I_Pssub,ijk) & ! [prod] s->v
           + w(I_Pgsub,ijk) & ! [prod] g->v
           - w(I_Psdep,ijk) & ! [loss] v->s
           - w(I_Pgdep,ijk)   ! [loss] v->g

       fac_sw = 0.5_RP + sign( 0.5_RP, net+EPS ) ! if production > loss , fac_sw=1
       fac    = (          fac_sw ) &
              + ( 1.0_RP - fac_sw ) * min( -w(I_dqv_dt,ijk)/(net-fac_sw), 1.0_RP ) ! loss limiter

       fac    = (          w(I_delta3,ijk) ) * fac &
              + ( 1.0_RP - w(I_delta3,ijk) )

       w(I_Prevp,ijk) = w(I_Prevp,ijk) * fac
       w(I_Pssub,ijk) = w(I_Pssub,ijk) * fac
       w(I_Pgsub,ijk) = w(I_Pgsub,ijk) * fac
       w(I_Psdep,ijk) = w(I_Psdep,ijk) * fac
       w(I_Pgdep,ijk) = w(I_Pgdep,ijk) * fac

       ! [QS]
       net = w(I_Psdep  ,ijk) & ! [prod] v->s
           + w(I_Psacw  ,ijk) & ! [prod] c->s
           + w(I_Psfw   ,ijk) & ! [prod] c->s
           + w(I_Piacr_s,ijk) & ! [prod] r->s
           + w(I_Psacr_s,ijk) & ! [prod] r->s
           + w(I_Psaut  ,ijk) & ! [prod] i->s
           + w(I_Praci_s,ijk) & ! [prod] i->s
           + w(I_Psaci  ,ijk) & ! [prod] i->s
           + w(I_Psfi   ,ijk) & ! [prod] i->s
           - w(I_Pssub  ,ijk) & ! [loss] s->v
           - w(I_Pgaut  ,ijk) & ! [loss] s->g
           - w(I_Pracs  ,ijk) & ! [loss] s->g
           - w(I_Pgacs  ,ijk)   ! [loss] s->g

       fac_sw = 0.5_RP + sign( 0.5_RP, net+EPS ) ! if production > loss , fac_sw=1
       fac    = (          fac_sw ) &
              + ( 1.0_RP - fac_sw ) * min( -w(I_dqs_dt,ijk)/(net-fac_sw), 1.0_RP ) ! loss limiter

       w(I_Psdep  ,ijk) = w(I_Psdep  ,ijk) * fac
       w(I_Psacw  ,ijk) = w(I_Psacw  ,ijk) * fac
       w(I_Psfw   ,ijk) = w(I_Psfw   ,ijk) * fac
       w(I_Piacr_s,ijk) = w(I_Piacr_s,ijk) * fac
       w(I_Psacr_s,ijk) = w(I_Psacr_s,ijk) * fac
       w(I_Psaut  ,ijk) = w(I_Psaut  ,ijk) * fac
       w(I_Praci_s,ijk) = w(I_Praci_s,ijk) * fac
       w(I_Psaci  ,ijk) = w(I_Psaci  ,ijk) * fac
       w(I_Psfi   ,ijk) = w(I_Psfi   ,ijk) * fac
       w(I_Pssub  ,ijk) = w(I_Pssub  ,ijk) * fac
       w(I_Pracs  ,ijk) = w(I_Pracs  ,ijk) * fac
       w(I_Pgaut  ,ijk) = w(I_Pgaut  ,ijk) * fac
       w(I_Pgacs  ,ijk) = w(I_Pgacs  ,ijk) * fac

       ! [QG]
       net = w(I_Pgdep  ,ijk) & ! [prod] v->g
           + w(I_Pgacw  ,ijk) & ! [prod] c->g
           + w(I_Piacr_g,ijk) & ! [prod] r->g
           + w(I_Psacr_g,ijk) & ! [prod] r->g
           + w(I_Pgacr  ,ijk) & ! [prod] r->g
           + w(I_Pgfrz  ,ijk) & ! [prod] r->g
           + w(I_Praci_g,ijk) & ! [prod] i->g
           + w(I_Pgaci  ,ijk) & ! [prod] i->g
           + w(I_Pgaut  ,ijk) & ! [prod] s->g
           + w(I_Pracs  ,ijk) & ! [prod] s->g
           + w(I_Pgacs  ,ijk) & ! [prod] s->g
           - w(I_Pgsub  ,ijk)   ! [loss] g->v

       fac_sw = 0.5_RP + sign( 0.5_RP, net+EPS ) ! if production > loss , fac_sw=1
       fac    = (          fac_sw ) &
              + ( 1.0_RP - fac_sw ) * min( -w(I_dqg_dt,ijk)/(net-fac_sw), 1.0_RP ) ! loss limiter

       w(I_Pgdep  ,ijk) = w(I_Pgdep  ,ijk) * fac
       w(I_Pgacw  ,ijk) = w(I_Pgacw  ,ijk) * fac
       w(I_Piacr_g,ijk) = w(I_Piacr_g,ijk) * fac
       w(I_Psacr_g,ijk) = w(I_Psacr_g,ijk) * fac
       w(I_Pgacr  ,ijk) = w(I_Pgacr  ,ijk) * fac
       w(I_Pgfrz  ,ijk) = w(I_Pgfrz  ,ijk) * fac
       w(I_Praci_g,ijk) = w(I_Praci_g,ijk) * fac
       w(I_Pgaci  ,ijk) = w(I_Pgaci  ,ijk) * fac
       w(I_Pgaut  ,ijk) = w(I_Pgaut  ,ijk) * fac
       w(I_Pracs  ,ijk) = w(I_Pracs  ,ijk) * fac
       w(I_Pgacs  ,ijk) = w(I_Pgacs  ,ijk) * fac
       w(I_Pgsub  ,ijk) = w(I_Pgsub  ,ijk) * fac
    enddo



    do indirect = 1, ijk_cold
       ijk = index_cold(indirect)

       ! re-calc net production & loss
       tend(I_QC,ijk) = &                  ! [prod]
                      - w(I_Praut  ,ijk) & ! [loss] c->r
                      - w(I_Pracw  ,ijk) & ! [loss] c->r
                      - w(I_Psacw  ,ijk) & ! [loss] c->s
                      - w(I_Psfw   ,ijk) & ! [loss] c->s
                      - w(I_Pgacw  ,ijk)   ! [loss] c->g
       tend(I_QR,ijk) = w(I_Praut  ,ijk) & ! [prod] c->r
                      + w(I_Pracw  ,ijk) & ! [prod] c->r
                      - w(I_Prevp  ,ijk) & ! [loss] r->v
                      - w(I_Piacr_s,ijk) & ! [loss] r->s
                      - w(I_Psacr_s,ijk) & ! [loss] r->s
                      - w(I_Piacr_g,ijk) & ! [loss] r->g
                      - w(I_Psacr_g,ijk) & ! [loss] r->g
                      - w(I_Pgacr  ,ijk) & ! [loss] r->g
                      - w(I_Pgfrz  ,ijk)   ! [loss] r->g
       tend(I_QI,ijk) = &                  ! [prod]
                      - w(I_Psaut  ,ijk) & ! [loss] i->s
                      - w(I_Praci_s,ijk) & ! [loss] i->s
                      - w(I_Psaci  ,ijk) & ! [loss] i->s
                      - w(I_Psfi   ,ijk) & ! [loss] i->s
                      - w(I_Praci_g,ijk) & ! [loss] i->g
                      - w(I_Pgaci  ,ijk)   ! [loss] i->g
       tend(I_QS,ijk) = w(I_Psdep  ,ijk) & ! [prod] v->s
                      + w(I_Psacw  ,ijk) & ! [prod] c->s
                      + w(I_Psfw   ,ijk) & ! [prod] c->s
                      + w(I_Piacr_s,ijk) & ! [prod] r->s
                      + w(I_Psacr_s,ijk) & ! [prod] r->s
                      + w(I_Psaut  ,ijk) & ! [prod] i->s
                      + w(I_Psaci  ,ijk) & ! [prod] i->s
                      + w(I_Praci_s,ijk) & ! [prod] i->s
                      + w(I_Psfi   ,ijk) & ! [prod] i->s
                      - w(I_Pssub  ,ijk) & ! [loss] s->v
                      - w(I_Pgaut  ,ijk) & ! [loss] s->g
                      - w(I_Pracs  ,ijk) & ! [loss] s->g
                      - w(I_Pgacs  ,ijk)   ! [loss] s->g
       tend(I_QG,ijk) = w(I_Pgdep  ,ijk) & ! [prod] v->g
                      + w(I_Pgacw  ,ijk) & ! [prod] c->g
                      + w(I_Piacr_g,ijk) & ! [prod] r->g
                      + w(I_Psacr_g,ijk) & ! [prod] r->g
                      + w(I_Pgacr  ,ijk) & ! [prod] r->g
                      + w(I_Pgfrz  ,ijk) & ! [prod] r->g
                      + w(I_Praci_g,ijk) & ! [prod] i->g
                      + w(I_Pgaci  ,ijk) & ! [prod] i->g
                      + w(I_Pgaut  ,ijk) & ! [prod] s->g
                      + w(I_Pracs  ,ijk) & ! [prod] s->g
                      + w(I_Pgacs  ,ijk) & ! [prod] s->g
                      - w(I_Pgsub  ,ijk)   ! [loss] g->v

       tend(I_QC,ijk) = max( tend(I_QC,ijk), -w(I_dqc_dt,ijk) )
       tend(I_QR,ijk) = max( tend(I_QR,ijk), -w(I_dqr_dt,ijk) )
       tend(I_QI,ijk) = max( tend(I_QI,ijk), -w(I_dqi_dt,ijk) )
       tend(I_QS,ijk) = max( tend(I_QS,ijk), -w(I_dqs_dt,ijk) )
       tend(I_QG,ijk) = max( tend(I_QG,ijk), -w(I_dqg_dt,ijk) )

       tend(I_QV,ijk) = - ( tend(I_QC,ijk) &
                          + tend(I_QR,ijk) &
                          + tend(I_QI,ijk) &
                          + tend(I_QS,ijk) &
                          + tend(I_QG,ijk) )
    enddo

    ! mass & energy update
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       ijk = ( j - JS ) * KMAX * IMAX &
           + ( i - IS ) * KMAX        &
           + ( k - KS )               &
           + 1

       QTRC_t(k,i,j,I_QV) = QTRC_t(k,i,j,I_QV) + tend(I_QV,ijk)
       QTRC_t(k,i,j,I_QC) = QTRC_t(k,i,j,I_QC) + tend(I_QC,ijk)
       QTRC_t(k,i,j,I_QR) = QTRC_t(k,i,j,I_QR) + tend(I_QR,ijk)
       QTRC_t(k,i,j,I_QI) = QTRC_t(k,i,j,I_QI) + tend(I_QI,ijk)
       QTRC_t(k,i,j,I_QS) = QTRC_t(k,i,j,I_QS) + tend(I_QS,ijk)
       QTRC_t(k,i,j,I_QG) = QTRC_t(k,i,j,I_QG) + tend(I_QG,ijk)

       RHOE_t(k,i,j) = RHOE_t(k,i,j) - DENS0(k,i,j) * ( LHV00 * tend(I_QV,ijk) &
                                                      - LHF00 * tend(I_QI,ijk) &
                                                      - LHF00 * tend(I_QS,ijk) &
                                                      - LHF00 * tend(I_QG,ijk) )

       QTRC0(k,i,j,I_QV) = QTRC0(k,i,j,I_QV) + QTRC_t(k,i,j,I_QV) * dt
       QTRC0(k,i,j,I_QC) = QTRC0(k,i,j,I_QC) + QTRC_t(k,i,j,I_QC) * dt
       QTRC0(k,i,j,I_QR) = QTRC0(k,i,j,I_QR) + QTRC_t(k,i,j,I_QR) * dt
       QTRC0(k,i,j,I_QI) = QTRC0(k,i,j,I_QI) + QTRC_t(k,i,j,I_QI) * dt
       QTRC0(k,i,j,I_QS) = QTRC0(k,i,j,I_QS) + QTRC_t(k,i,j,I_QS) * dt
       QTRC0(k,i,j,I_QG) = QTRC0(k,i,j,I_QG) + QTRC_t(k,i,j,I_QG) * dt

       RHOE0(k,i,j) = RHOE0(k,i,j) + RHOE_t(k,i,j) * dt
    enddo
    enddo
    enddo

    do ip = 1, w_nmax
       if ( w_histid(ip) > 0 ) then
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             ijk = ( j - JS ) * KMAX * IMAX &
                 + ( i - IS ) * KMAX        &
                 + ( k - KS )               &
                 + 1

             work3D(k,i,j) = w(ip,ijk)
          enddo
          enddo
          enddo

!          if( IO_L ) write(IO_FID_LOG,*) w_name(ip), "MAX/MIN:", &
!                     maxval(work3D(KS:KE,IS:IE,JS:JE)), minval(work3D(KS:KE,IS:IE,JS:JE))

          call HIST_put( w_histid(ip), work3D(:,:,:), dt, w_zinterp(ip) )
       endif
    enddo

!    if( IO_L ) write(IO_FID_LOG,*) "tend QV MAX/MIN:", maxval(tend(I_QV,:)), minval(tend(I_QV,:))
!    if( IO_L ) write(IO_FID_LOG,*) "tend QC MAX/MIN:", maxval(tend(I_QC,:)), minval(tend(I_QC,:))
!    if( IO_L ) write(IO_FID_LOG,*) "tend QR MAX/MIN:", maxval(tend(I_QR,:)), minval(tend(I_QR,:))
!    if( IO_L ) write(IO_FID_LOG,*) "tend QI MAX/MIN:", maxval(tend(I_QI,:)), minval(tend(I_QI,:))
!    if( IO_L ) write(IO_FID_LOG,*) "tend QS MAX/MIN:", maxval(tend(I_QS,:)), minval(tend(I_QS,:))
!    if( IO_L ) write(IO_FID_LOG,*) "tend QG MAX/MIN:", maxval(tend(I_QG,:)), minval(tend(I_QG,:))
!
!    if( IO_L ) write(IO_FID_LOG,*) "QV MAX/MIN:", maxval(QTRC0(:,:,:,I_QV)), minval(QTRC0(:,:,:,I_QV))
!    if( IO_L ) write(IO_FID_LOG,*) "QC MAX/MIN:", maxval(QTRC0(:,:,:,I_QC)), minval(QTRC0(:,:,:,I_QC))
!    if( IO_L ) write(IO_FID_LOG,*) "QR MAX/MIN:", maxval(QTRC0(:,:,:,I_QR)), minval(QTRC0(:,:,:,I_QR))
!    if( IO_L ) write(IO_FID_LOG,*) "QI MAX/MIN:", maxval(QTRC0(:,:,:,I_QI)), minval(QTRC0(:,:,:,I_QI))
!    if( IO_L ) write(IO_FID_LOG,*) "QS MAX/MIN:", maxval(QTRC0(:,:,:,I_QS)), minval(QTRC0(:,:,:,I_QS))
!    if( IO_L ) write(IO_FID_LOG,*) "QG MAX/MIN:", maxval(QTRC0(:,:,:,I_QG)), minval(QTRC0(:,:,:,I_QG))

!    if( IO_L ) write(IO_FID_LOG,*) "QV_t MAX/MIN:", maxval(QTRC_t(:,:,:,I_QV)), minval(QTRC_t(:,:,:,I_QV))
!    if( IO_L ) write(IO_FID_LOG,*) "QC_t MAX/MIN:", maxval(QTRC_t(:,:,:,I_QC)), minval(QTRC_t(:,:,:,I_QC))
!    if( IO_L ) write(IO_FID_LOG,*) "QR_t MAX/MIN:", maxval(QTRC_t(:,:,:,I_QR)), minval(QTRC_t(:,:,:,I_QR))
!    if( IO_L ) write(IO_FID_LOG,*) "QI_t MAX/MIN:", maxval(QTRC_t(:,:,:,I_QI)), minval(QTRC_t(:,:,:,I_QI))
!    if( IO_L ) write(IO_FID_LOG,*) "QS_t MAX/MIN:", maxval(QTRC_t(:,:,:,I_QS)), minval(QTRC_t(:,:,:,I_QS))
!    if( IO_L ) write(IO_FID_LOG,*) "QG_t MAX/MIN:", maxval(QTRC_t(:,:,:,I_QG)), minval(QTRC_t(:,:,:,I_QG))

!    if( IO_L ) write(IO_FID_LOG,*) "RHOE0  MAX/MIN:", maxval(RHOE0(:,:,:)), minval(RHOE0(:,:,:))
!    if( IO_L ) write(IO_FID_LOG,*) "RHOE_t MAX/MIN:", maxval(RHOE_t(:,:,:)), minval(RHOE_t(:,:,:))

    call PROF_rapend  ('MP_tomita08')

    return
  end subroutine MP_tomita08

  !-----------------------------------------------------------------------------
  !> Lin-type cold rain microphysics (terminal velocity)
  !-----------------------------------------------------------------------------
  subroutine MP_tomita08_vterm( &
       vterm, &
       DENS0,  &
       QTRC0   )
    implicit none

    real(RP), intent(inout) :: vterm(KA,IA,JA,QA)
    real(RP), intent(in)    :: DENS0(KA,IA,JA)
    real(RP), intent(in)    :: QTRC0(KA,IA,JA,QA)

    real(RP) :: dens
    real(RP) :: q(QA)

    real(RP) :: rho_fact ! density factor

    real(RP) :: RLMDr, RLMDs, RLMDg
    real(RP) :: RLMDr_dr, RLMDs_ds, RLMDg_dg

    real(RP) :: zerosw
    integer  :: k, i, j, iq
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       ! store to work
       dens = DENS0(k,i,j)
       do iq = I_QV, I_QG
          q(iq) = QTRC0(k,i,j,iq)
       enddo

       rho_fact = sqrt( dens00 / dens )

       ! slope parameter lambda
       zerosw = 0.5_RP - sign(0.5_RP, q(I_QR) - 1.E-12_RP )
       RLMDr = sqrt(sqrt( dens * q(I_QR) / ( Ar * N0r * GAM_1br ) + zerosw )) * ( 1.0_RP - zerosw )

       zerosw = 0.5_RP - sign(0.5_RP, q(I_QS) - 1.E-12_RP )
       RLMDs = sqrt(sqrt( dens * q(I_QS) / ( As * N0s * GAM_1bs ) + zerosw )) * ( 1.0_RP - zerosw )

       zerosw = 0.5_RP - sign(0.5_RP, q(I_QG) - 1.E-12_RP )
       RLMDg = sqrt(sqrt( dens * q(I_QG) / ( Ag * N0g * GAM_1bg ) + zerosw )) * ( 1.0_RP - zerosw )

       RLMDr_dr = sqrt( RLMDr )       ! **Dr
       RLMDs_ds = sqrt( sqrt(RLMDs) ) ! **Ds
       RLMDg_dg = sqrt( RLMDg )       ! **Dg

       !---< terminal velocity >
       zerosw = 0.5_RP + sign(0.5_RP, q(I_QI) - 1.E-8_RP )
       vterm(k,i,j,I_QI) = -3.29_RP * ( dens * q(I_QI) * zerosw )**0.16_RP
       vterm(k,i,j,I_QR) = -Cr * rho_fact * GAM_1brdr / GAM_1br * RLMDr_dr
       vterm(k,i,j,I_QS) = -Cs * rho_fact * GAM_1bsds / GAM_1bs * RLMDs_ds
       vterm(k,i,j,I_QG) = -Cg * rho_fact * GAM_1bgdg / GAM_1bg * RLMDg_dg
    enddo
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
  end subroutine ATMOS_PHY_MP_tomita08_CloudFraction

  !-----------------------------------------------------------------------------
  !> Calculate Effective Radius
  subroutine ATMOS_PHY_MP_tomita08_EffectiveRadius( &
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

    real(RP) :: dens
    real(RP) :: q(QA)
    real(RP) :: RLMDr, RLMDs, RLMDg

    real(RP) :: zerosw
    integer  :: k, i, j, iq
    !---------------------------------------------------------------------------

    Re(:,:,:,I_mp_QC) =   8.E-6_RP
    Re(:,:,:,I_mp_QI) =  20.E-6_RP

    ! Effective radius is defined by r3m/r2m=1.5/lambda.
    do j  = JS, JE
    do i  = IS, IE
    do k  = KS, KE
       ! store to work
       dens = DENS0(k,i,j)
       do iq = I_QV, I_QG
          q(iq) = QTRC0(k,i,j,iq)
       enddo

       ! slope parameter lambda
       zerosw = 0.5_RP - sign(0.5_RP, q(I_QR) - 1.E-12_RP )
       RLMDr = sqrt(sqrt( dens * q(I_QR) / ( Ar * N0r * GAM_1br ) + zerosw )) * ( 1.0_RP - zerosw )

       zerosw = 0.5_RP - sign(0.5_RP, q(I_QS) - 1.E-12_RP )
       RLMDs = sqrt(sqrt( dens * q(I_QS) / ( As * N0s * GAM_1bs ) + zerosw )) * ( 1.0_RP - zerosw )

       zerosw = 0.5_RP - sign(0.5_RP, q(I_QG) - 1.E-12_RP )
       RLMDg = sqrt(sqrt( dens * q(I_QG) / ( Ag * N0g * GAM_1bg ) + zerosw )) * ( 1.0_RP - zerosw )

       Re(k,i,j,I_mp_QR) = RLMDr * 1.5_RP
       Re(k,i,j,I_mp_QS) = RLMDs * 1.5_RP
       Re(k,i,j,I_mp_QG) = RLMDg * 1.5_RP
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_MP_tomita08_EffectiveRadius
  !-----------------------------------------------------------------------------
  !> Calculate mixing ratio of each category
  subroutine ATMOS_PHY_MP_tomita08_Mixingratio( &
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

  end subroutine ATMOS_PHY_MP_tomita08_Mixingratio
  !-----------------------------------------------------------------------------
end module scale_atmos_phy_mp_tomita08
!-------------------------------------------------------------------------------
