!-------------------------------------------------------------------------------
!> module ATMOSPHERIC Variables
!!
!! @par Description
!!          Container for atmospheric variables
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2011-11-11 (H.Yashiro)   [new]
!! @li      2012-03-23 (H.Yashiro)   [mod] Explicit index parameter inclusion
!! @li      2012-06-13 (S.Nishizawa) [mod] follows the  change of mod_hist
!!
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
#include "macro_thermodyn.h"
module mod_atmos_vars
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_debug
  use scale_grid_index
  use scale_index
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_vars_setup
  public :: ATMOS_vars_fillhalo
  public :: ATMOS_vars_restart_read
  public :: ATMOS_vars_restart_write
  public :: ATMOS_vars_restart_check
  public :: ATMOS_vars_history_setpres
  public :: ATMOS_vars_history
  public :: ATMOS_vars_total
  public :: ATMOS_vars_diagnostics
  public :: ATMOS_vars_monitor

  public :: ATMOS_vars_restart_create
  public :: ATMOS_vars_restart_open
  public :: ATMOS_vars_restart_def_var
  public :: ATMOS_vars_restart_enddef
  public :: ATMOS_vars_restart_close

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical,               public :: ATMOS_RESTART_OUTPUT                = .false.         !< Output restart file?

  character(len=H_LONG), public :: ATMOS_RESTART_IN_BASENAME           = ''              !< Basename of the input  file
  logical,               public :: ATMOS_RESTART_IN_POSTFIX_TIMELABEL  = .false.         !< Add timelabel to the basename of input  file?
  character(len=H_LONG), public :: ATMOS_RESTART_OUT_BASENAME          = ''              !< Basename of the output file
  logical,               public :: ATMOS_RESTART_OUT_POSTFIX_TIMELABEL = .true.          !< Add timelabel to the basename of output file?
  character(len=H_MID),  public :: ATMOS_RESTART_OUT_TITLE             = 'ATMOS restart' !< Title    of the output file
  character(len=H_MID),  public :: ATMOS_RESTART_OUT_DTYPE             = 'DEFAULT'       !< REAL4 or REAL8

  logical,               public :: ATMOS_RESTART_CHECK                 = .false.         !< Check value consistency?
  character(len=H_LONG), public :: ATMOS_RESTART_CHECK_BASENAME        = 'restart_check'
  real(RP),              public :: ATMOS_RESTART_CHECK_CRITERION       = 1.E-6_RP

  ! prognostic variables
  real(RP), public, target, allocatable :: DENS(:,:,:)   ! Density     [kg/m3]
  real(RP), public, target, allocatable :: MOMZ(:,:,:)   ! momentum z  [kg/m2/s]
  real(RP), public, target, allocatable :: MOMX(:,:,:)   ! momentum x  [kg/m2/s]
  real(RP), public, target, allocatable :: MOMY(:,:,:)   ! momentum y  [kg/m2/s]
  real(RP), public, target, allocatable :: RHOT(:,:,:)   ! DENS * POTT [K*kg/m3]
  real(RP), public, target, allocatable :: QTRC(:,:,:,:) ! ratio of mass of tracer to total mass[kg/kg]

  real(RP), public, target, allocatable :: DENS_avw(:,:,:)
  real(RP), public, target, allocatable :: MOMZ_avw(:,:,:)
  real(RP), public, target, allocatable :: MOMX_avw(:,:,:)
  real(RP), public, target, allocatable :: MOMY_avw(:,:,:)
  real(RP), public, target, allocatable :: RHOT_avw(:,:,:)
  real(RP), public, target, allocatable :: QTRC_avw(:,:,:,:)

  real(RP), public, pointer             :: DENS_av(:,:,:)
  real(RP), public, pointer             :: MOMZ_av(:,:,:)
  real(RP), public, pointer             :: MOMX_av(:,:,:)
  real(RP), public, pointer             :: MOMY_av(:,:,:)
  real(RP), public, pointer             :: RHOT_av(:,:,:)
  real(RP), public, pointer             :: QTRC_av(:,:,:,:)

  ! tendency by physical processes
  real(RP), public, allocatable :: DENS_tp(:,:,:)
  real(RP), public, allocatable :: MOMZ_tp(:,:,:)
  real(RP), public, allocatable :: MOMX_tp(:,:,:)
  real(RP), public, allocatable :: MOMY_tp(:,:,:)
  real(RP), public, allocatable :: RHOT_tp(:,:,:)
  real(RP), public, allocatable :: RHOQ_tp(:,:,:,:)

  ! diagnostic variables
  real(RP), public, allocatable :: TEMP(:,:,:)   ! temperature [K]
  real(RP), public, allocatable :: PRES(:,:,:)   ! pressure    [Pa=J/m3]
  real(RP), public, allocatable :: W   (:,:,:)   ! velocity w  [m/s]
  real(RP), public, allocatable :: U   (:,:,:)   ! velocity u  [m/s]
  real(RP), public, allocatable :: V   (:,:,:)   ! velocity v  [m/s]
  real(RP), public, allocatable :: POTT(:,:,:)   ! potential temperature [K]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical,                private              :: ATMOS_VARS_CHECKRANGE = .false.
  real(RP),               private              :: ATMOS_VARS_CHECKCFL   = 0.0_RP

  integer,                private, parameter   :: VMAX   = 5       !< number of the variables
  character(len=H_SHORT), private              :: VAR_NAME(VMAX)
  character(len=H_MID),   private              :: VAR_DESC(VMAX)
  character(len=H_SHORT), private              :: VAR_UNIT(VMAX)
  integer,                private, allocatable :: VAR_ID(:)
  integer,                private              :: restart_fid = -1  ! file ID

  data VAR_NAME / 'DENS', &
                  'MOMZ', &
                  'MOMX', &
                  'MOMY', &
                  'RHOT'  /
  data VAR_DESC / 'density',    &
                  'momentum z', &
                  'momentum x', &
                  'momentum y', &
                  'rho * theta' /
  data VAR_UNIT / 'kg/m3',   &
                  'kg/m2/s', &
                  'kg/m2/s', &
                  'kg/m2/s', &
                  'kg/m3*K'  /

  ! history output of prognostic variables
  integer, private              :: VAR_HIST_id(VMAX)
  integer, private, allocatable :: AQ_HIST_id(:)

  ! history & monitor output of diagnostic variables
  integer, private, parameter :: AD_nmax = 70        ! number of diagnostic variables for history output

  integer, private, parameter :: I_W            =  1 ! velocity w at cell center
  integer, private, parameter :: I_U            =  2 ! velocity u at cell center
  integer, private, parameter :: I_V            =  3 ! velocity v at cell center
  integer, private, parameter :: I_POTT         =  4 ! potential temperature

  integer, private, parameter :: I_QDRY         =  5 ! ratio of dry air            to total mass
  integer, private, parameter :: I_QTOT         =  6 ! ratio of total tracer       to total mass
  integer, private, parameter :: I_QHYD         =  7 ! ratio of total hydrometeor  to total mass
  integer, private, parameter :: I_QLIQ         =  8 ! ratio of total liquid water to total mass
  integer, private, parameter :: I_QICE         =  9 ! ratio of total ice    water to total mass

  integer, private, parameter :: I_LWP          = 10 ! liquid water path
  integer, private, parameter :: I_IWP          = 11 ! ice water path
  integer, private, parameter :: I_PW           = 12 ! ice water path

  integer, private, parameter :: I_RTOT         = 13 ! total gas constant
  integer, private, parameter :: I_CPTOT        = 14 ! total heat capacity (constant pressure)
  integer, private, parameter :: I_PRES         = 15 ! pressure
  integer, private, parameter :: I_TEMP         = 16 ! temperature

  integer, private, parameter :: I_POTL         = 17 ! liquid water potential temperature
  integer, private, parameter :: I_RHA          = 18 ! relative humidity (liquid+ice)
  integer, private, parameter :: I_RHL          = 19 ! relative humidity against to liquid
  integer, private, parameter :: I_RHI          = 20 ! relative humidity against to ice

  integer, private, parameter :: I_VOR          = 21 ! vertical vorticity
  integer, private, parameter :: I_DIV          = 22 ! divergence
  integer, private, parameter :: I_HDIV         = 23 ! horizontal divergence

  integer, private, parameter :: I_DENS_PRIM    = 24 ! prime term of density
  integer, private, parameter :: I_W_PRIM       = 25 ! prime term of w
  integer, private, parameter :: I_U_PRIM       = 26 ! prime term of u
  integer, private, parameter :: I_V_PRIM       = 27 ! prime term of v
  integer, private, parameter :: I_POTT_PRIM    = 28 ! prime term of potential temperature
  integer, private, parameter :: I_W_PRIM2      = 29 ! variance of w
  integer, private, parameter :: I_PT_W_PRIM    = 30 ! resolved scale heat flux
  integer, private, parameter :: I_W_PRIM3      = 31 ! skewness of w
  integer, private, parameter :: I_TKE_RS       = 32 ! resolved scale TKE

  integer, private, parameter :: I_ENGP         = 33 ! potential energy
  integer, private, parameter :: I_ENGK         = 34 ! kinetic   energy
  integer, private, parameter :: I_ENGI         = 35 ! internal  energy
  integer, private, parameter :: I_ENGT         = 36 ! total     energy

  integer, private, parameter :: I_ENGSFC_SH    = 37
  integer, private, parameter :: I_ENGSFC_LH    = 38
  integer, private, parameter :: I_ENGSFC_RD    = 39
  integer, private, parameter :: I_ENGTOA_RD    = 40

  integer, private, parameter :: I_ENGSFC_LW_up = 41
  integer, private, parameter :: I_ENGSFC_LW_dn = 42
  integer, private, parameter :: I_ENGSFC_SW_up = 43
  integer, private, parameter :: I_ENGSFC_SW_dn = 44

  integer, private, parameter :: I_ENGTOA_LW_up = 45
  integer, private, parameter :: I_ENGTOA_LW_dn = 46
  integer, private, parameter :: I_ENGTOA_SW_up = 47
  integer, private, parameter :: I_ENGTOA_SW_dn = 48

  integer, private, parameter :: I_ENGFLXT      = 49

  integer, private, parameter :: I_EVAP         = 50
  integer, private, parameter :: I_PRCP         = 51

  integer, private, parameter :: I_DENS_MEAN    = 52
  integer, private, parameter :: I_W_MEAN       = 53
  integer, private, parameter :: I_U_MEAN       = 54
  integer, private, parameter :: I_V_MEAN       = 55
  integer, private, parameter :: I_POTT_MEAN    = 56
  integer, private, parameter :: I_T_MEAN       = 57

  integer, private, parameter :: I_QV_MEAN      = 58
  integer, private, parameter :: I_QHYD_MEAN    = 59
  integer, private, parameter :: I_QLIQ_MEAN    = 60
  integer, private, parameter :: I_QICE_MEAN    = 61

  integer, private, parameter :: I_QSAT         = 62

  integer, private, parameter :: I_Uabs         = 63

  integer, private, parameter :: I_CAPE         = 64
  integer, private, parameter :: I_CIN          = 65
  integer, private, parameter :: I_LCL          = 66
  integer, private, parameter :: I_LFC          = 67
  integer, private, parameter :: I_LNB          = 68

  integer, private, parameter :: I_PBLH         = 69
  integer, private, parameter :: I_MSE          = 70

  integer, private            :: AD_HIST_id (AD_nmax)
  integer, private            :: AD_PREP_sw (AD_nmax)
  integer, private            :: AD_MONIT_id(AD_nmax)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_vars_setup
    use scale_process, only: &
       PRC_MPIstop
    use scale_history, only: &
       HIST_reg
    use scale_monitor, only: &
       MONIT_reg
    use mod_atmos_admin, only: &
       ATMOS_USE_AVERAGE
    use mod_atmos_dyn_vars, only: &
       ATMOS_DYN_vars_setup
    use mod_atmos_phy_mp_vars, only: &
       ATMOS_PHY_MP_vars_setup
    use mod_atmos_phy_ae_vars, only: &
       ATMOS_PHY_AE_vars_setup
    use mod_atmos_phy_ch_vars, only: &
       ATMOS_PHY_CH_vars_setup
    use mod_atmos_phy_rd_vars, only: &
       ATMOS_PHY_RD_vars_setup
    use mod_atmos_phy_sf_vars, only: &
       ATMOS_PHY_SF_vars_setup
    use mod_atmos_phy_tb_vars, only: &
       ATMOS_PHY_TB_vars_setup
    use mod_atmos_phy_cp_vars, only: &
       ATMOS_PHY_CP_vars_setup
    implicit none

    NAMELIST / PARAM_ATMOS_VARS / &
       ATMOS_RESTART_IN_BASENAME,           &
       ATMOS_RESTART_IN_POSTFIX_TIMELABEL,  &
       ATMOS_RESTART_OUTPUT,                &
       ATMOS_RESTART_OUT_BASENAME,          &
       ATMOS_RESTART_OUT_POSTFIX_TIMELABEL, &
       ATMOS_RESTART_OUT_TITLE,             &
       ATMOS_RESTART_OUT_DTYPE,             &
       ATMOS_RESTART_CHECK,                 &
       ATMOS_RESTART_CHECK_BASENAME,        &
       ATMOS_RESTART_CHECK_CRITERION,       &
       ATMOS_VARS_CHECKRANGE,               &
       ATMOS_VARS_CHECKCFL

    integer :: ierr
    integer :: iv, iq
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[VARS] / Categ[ATMOS] / Origin[SCALE-RM]'

    allocate( DENS(KA,IA,JA)    )
    allocate( MOMZ(KA,IA,JA)    )
    allocate( MOMX(KA,IA,JA)    )
    allocate( MOMY(KA,IA,JA)    )
    allocate( RHOT(KA,IA,JA)    )
    allocate( QTRC(KA,IA,JA,max(QA,1)) )

    if ( ATMOS_USE_AVERAGE ) then
       allocate( DENS_avw(KA,IA,JA)    )
       allocate( MOMZ_avw(KA,IA,JA)    )
       allocate( MOMX_avw(KA,IA,JA)    )
       allocate( MOMY_avw(KA,IA,JA)    )
       allocate( RHOT_avw(KA,IA,JA)    )
       allocate( QTRC_avw(KA,IA,JA,max(QA,1)) )

       DENS_av => DENS_avw
       MOMZ_av => MOMZ_avw
       MOMX_av => MOMX_avw
       MOMY_av => MOMY_avw
       RHOT_av => RHOT_avw
       QTRC_av => QTRC_avw
    else
       DENS_av => DENS
       MOMZ_av => MOMZ
       MOMX_av => MOMX
       MOMY_av => MOMY
       RHOT_av => RHOT
       QTRC_av => QTRC
    endif

    allocate( DENS_tp(KA,IA,JA)    )
    allocate( MOMZ_tp(KA,IA,JA)    )
    allocate( MOMX_tp(KA,IA,JA)    )
    allocate( MOMY_tp(KA,IA,JA)    )
    allocate( RHOT_tp(KA,IA,JA)    )
    allocate( RHOQ_tp(KA,IA,JA,max(QA,1)) )

    allocate( TEMP(KA,IA,JA) )
    allocate( PRES(KA,IA,JA) )
    allocate( W   (KA,IA,JA) )
    allocate( U   (KA,IA,JA) )
    allocate( V   (KA,IA,JA) )
    allocate( POTT(KA,IA,JA) )

    MOMZ(1:KS-1,:,:) = 0.0_RP
    MOMZ(KE:KA,:,:) = 0.0_RP

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_VARS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_VARS. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_ATMOS_VARS)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** List of prognostic variables (ATMOS) ***'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A15,A,A32,3(A))') &
               '***       |','VARNAME        ','|', 'DESCRIPTION                     ','[', 'UNIT            ',']'
    do iv = 1, VMAX
       if( IO_L ) write(IO_FID_LOG,'(1x,A,i3,A,A15,A,A32,3(A))') &
                  '*** NO.',iv,'|',VAR_NAME(iv),'|', VAR_DESC(iv),'[', VAR_UNIT(iv),']'
    enddo
    do iq = 1, QA
       if( IO_L ) write(IO_FID_LOG,'(1x,A,i3,A,A15,A,A32,3(A))') &
                  '*** NO.',5+iq,'|',TRACER_NAME(iq),'|', TRACER_DESC(iq),'[', TRACER_UNIT(iq),']'
    enddo

    if( IO_L ) write(IO_FID_LOG,*)
    if ( ATMOS_RESTART_IN_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart input?  : YES, file = ', trim(ATMOS_RESTART_IN_BASENAME)
       if( IO_L ) write(IO_FID_LOG,*) '*** Add timelabel?  : ', ATMOS_RESTART_IN_POSTFIX_TIMELABEL
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart input?  : NO'
    endif
    if (       ATMOS_RESTART_OUTPUT             &
         .AND. ATMOS_RESTART_OUT_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart output? : YES, file = ', trim(ATMOS_RESTART_OUT_BASENAME)
       if( IO_L ) write(IO_FID_LOG,*) '*** Add timelabel?  : ', ATMOS_RESTART_OUT_POSTFIX_TIMELABEL
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart output? : NO'
       ATMOS_RESTART_OUTPUT = .false.
    endif

    if ( ATMOS_RESTART_CHECK_BASENAME == '' ) then
       ATMOS_RESTART_CHECK = .false.
    endif

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Check restart consistency?      : ', ATMOS_RESTART_CHECK
    if( IO_L ) write(IO_FID_LOG,*) '*** Check value range of variables? : ', ATMOS_VARS_CHECKRANGE
    if ( ATMOS_VARS_CHECKCFL > 0.0_RP ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Check CFL condition?            : YES'
       if( IO_L ) write(IO_FID_LOG,*) '*** Limit of Courant number         : ', ATMOS_VARS_CHECKCFL
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Check CFL condition?            : NO'
    endif

    call ATMOS_DYN_vars_setup
    call ATMOS_PHY_MP_vars_setup
    call ATMOS_PHY_AE_vars_setup
    call ATMOS_PHY_CH_vars_setup
    call ATMOS_PHY_RD_vars_setup
    call ATMOS_PHY_SF_vars_setup
    call ATMOS_PHY_TB_vars_setup
    call ATMOS_PHY_CP_vars_setup



    !##### todo: the part below should be moved to the mod_atmos_diag #####

    allocate( AQ_HIST_id (max(QA,1)))

    VAR_HIST_id(:) = -1
    AQ_HIST_id (:) = -1
    AD_HIST_id (:) = -1
    AD_MONIT_id(:) = -1
    AD_PREP_sw (:) = -1

    call HIST_reg( VAR_HIST_id(I_DENS), VAR_NAME(I_DENS), VAR_DESC(I_DENS), VAR_UNIT(I_DENS), ndim=3 )
    call HIST_reg( VAR_HIST_id(I_MOMZ), VAR_NAME(I_MOMZ), VAR_DESC(I_MOMZ), VAR_UNIT(I_MOMZ), ndim=3, zdim='half' )
    call HIST_reg( VAR_HIST_id(I_MOMX), VAR_NAME(I_MOMX), VAR_DESC(I_MOMX), VAR_UNIT(I_MOMX), ndim=3, xdim='half' )
    call HIST_reg( VAR_HIST_id(I_MOMY), VAR_NAME(I_MOMY), VAR_DESC(I_MOMY), VAR_UNIT(I_MOMY), ndim=3, ydim='half' )
    call HIST_reg( VAR_HIST_id(I_RHOT), VAR_NAME(I_RHOT), VAR_DESC(I_RHOT), VAR_UNIT(I_RHOT), ndim=3 )
    do iq = 1, QA
       call HIST_reg( AQ_HIST_id(iq), TRACER_NAME(iq), TRACER_DESC(iq), TRACER_UNIT(iq), ndim=3 )
    enddo

    call HIST_reg( AD_HIST_id(I_W)        , 'W',         'velocity w',                     'm/s',    ndim=3 )
    call HIST_reg( AD_HIST_id(I_U)        , 'U',         'velocity u',                     'm/s',    ndim=3 )
    call HIST_reg( AD_HIST_id(I_V)        , 'V',         'velocity v',                     'm/s',    ndim=3 )
    call HIST_reg( AD_HIST_id(I_POTT)     , 'PT',        'potential temp.',                'K',      ndim=3 )

    call HIST_reg( AD_HIST_id(I_QDRY)     , 'QDRY',      'dry air',                        'kg/kg',  ndim=3 )
    call HIST_reg( AD_HIST_id(I_QTOT)     , 'QTOT',      'total water',                    'kg/kg',  ndim=3 )
    call HIST_reg( AD_HIST_id(I_QHYD)     , 'QHYD',      'total hydrometeors',             'kg/kg',  ndim=3 )
    call HIST_reg( AD_HIST_id(I_QLIQ)     , 'QLIQ',      'total liquid water',             'kg/kg',  ndim=3 )
    call HIST_reg( AD_HIST_id(I_QICE)     , 'QICE',      'total ice water',                'kg/kg',  ndim=3 )

    call HIST_reg( AD_HIST_id(I_LWP)      , 'LWP',       'liquid water path',              'g/m2',   ndim=2 )
    call HIST_reg( AD_HIST_id(I_IWP)      , 'IWP',       'ice water path',                 'g/m2',   ndim=2 )
    call HIST_reg( AD_HIST_id(I_PW )      , 'PW',        'precipitable water',             'g/m2',   ndim=2 )

    call HIST_reg( AD_HIST_id(I_RTOT)     , 'RTOT',      'Total gas constant',             'J/kg/K', ndim=3 )
    call HIST_reg( AD_HIST_id(I_CPTOT)    , 'CPTOT',     'Total heat capacity',            'J/kg/K', ndim=3 )
    call HIST_reg( AD_HIST_id(I_PRES)     , 'PRES',      'pressure',                       'Pa',     ndim=3 )
    call HIST_reg( AD_HIST_id(I_TEMP)     , 'T',         'temperature',                    'K',      ndim=3 )

    call HIST_reg( AD_HIST_id(I_POTL)     , 'LWPT',      'liq. potential temp.',           'K',      ndim=3 )
    call HIST_reg( AD_HIST_id(I_RHA)      , 'RHA',       'relative humidity(liq+ice)',     '%',      ndim=3 )
    call HIST_reg( AD_HIST_id(I_RHL)      , 'RH',        'relative humidity(liq)',         '%',      ndim=3 )
    call HIST_reg( AD_HIST_id(I_RHI)      , 'RHI',       'relative humidity(ice)',         '%',      ndim=3 )

    call HIST_reg( AD_HIST_id(I_VOR)      , 'VOR',       'vertical vorticity',             '1/s',    ndim=3 )
    call HIST_reg( AD_HIST_id(I_DIV)      , 'DIV',       'divergence',                     '1/s',    ndim=3 )
    call HIST_reg( AD_HIST_id(I_HDIV)     , 'HDIV',      'horizontal divergence',          '1/s',    ndim=3 )
    call HIST_reg( AD_HIST_id(I_Uabs)     , 'Uabs',      'absolute velocity',              'm/s',    ndim=3 )

    call HIST_reg( AD_HIST_id(I_CAPE)     , 'CAPE',      'convection avail. pot. energy',  'm2/s2',  ndim=2 )
    call HIST_reg( AD_HIST_id(I_CIN)      , 'CIN',       'convection inhibition',          'm2/s2',  ndim=2 )
    call HIST_reg( AD_HIST_id(I_LCL)      , 'LCL',       'lifted condensation level',      'm',      ndim=2 )
    call HIST_reg( AD_HIST_id(I_LFC)      , 'LFC',       'level of free convection',       'm',      ndim=2 )
    call HIST_reg( AD_HIST_id(I_LNB)      , 'LNB',       'level of neutral buoyancy',      'm',      ndim=2 )

    call HIST_reg( AD_HIST_id(I_PBLH)     , 'PBLH',      'PBL height',                     'm',      ndim=2 )
    call HIST_reg( AD_HIST_id(I_MSE)      , 'MSE',       'moist static energy',            'm2/s2',  ndim=3 )

    call HIST_reg( AD_HIST_id(I_DENS_MEAN), 'DENS_MEAN', 'horiz. mean of density',         'kg/m3',  ndim=1 )
    call HIST_reg( AD_HIST_id(I_W_MEAN)   , 'W_MEAN',    'horiz. mean of w',               'm/s',    ndim=1 )
    call HIST_reg( AD_HIST_id(I_U_MEAN)   , 'U_MEAN',    'horiz. mean of u',               'm/s',    ndim=1 )
    call HIST_reg( AD_HIST_id(I_V_MEAN)   , 'V_MEAN',    'horiz. mean of v',               'm/s',    ndim=1 )
    call HIST_reg( AD_HIST_id(I_POTT_MEAN), 'PT_MEAN',   'horiz. mean of pot.',            'K',      ndim=1 )
    call HIST_reg( AD_HIST_id(I_T_MEAN)   , 'T_MEAN',    'horiz. mean of t',               'K',      ndim=1 )
    call HIST_reg( AD_HIST_id(I_QV_MEAN)  , 'QV_MEAN',   'horiz. mean of QV',              '1',      ndim=1 )
    call HIST_reg( AD_HIST_id(I_QHYD_MEAN), 'QHYD_MEAN', 'horiz. mean of QHYD',            '1',      ndim=1 )
    call HIST_reg( AD_HIST_id(I_QLIQ_MEAN), 'QLIQ_MEAN', 'horiz. mean of QLIQ',            '1',      ndim=1 )
    call HIST_reg( AD_HIST_id(I_QICE_MEAN), 'QICE_MEAN', 'horiz. mean of QICE',            '1',      ndim=1 )

    call HIST_reg( AD_HIST_id(I_DENS_PRIM), 'DENS_PRIM', 'horiz. deviation of density',    'kg/m3',  ndim=3 )
    call HIST_reg( AD_HIST_id(I_W_PRIM   ), 'W_PRIM',    'horiz. deviation of w',          'm/s',    ndim=3 )
    call HIST_reg( AD_HIST_id(I_U_PRIM   ), 'U_PRIM',    'horiz. deviation of u',          'm/s',    ndim=3 )
    call HIST_reg( AD_HIST_id(I_V_PRIM   ), 'V_PRIM',    'horiz. deviation of v',          'm/s',    ndim=3 )
    call HIST_reg( AD_HIST_id(I_POTT_PRIM), 'PT_PRIM',   'horiz. deviation of pot. temp.', 'K',      ndim=3 )
    call HIST_reg( AD_HIST_id(I_W_PRIM2  ), 'W_PRIM2',   'variance of w',                  'm2/s2',  ndim=3 )
    call HIST_reg( AD_HIST_id(I_PT_W_PRIM), 'PT_W_PRIM', 'resolved scale heat flux',       'W/s',    ndim=3 )
    call HIST_reg( AD_HIST_id(I_W_PRIM3  ), 'W_PRIM3',   'skewness of w',                  'm3/s3',  ndim=3 )
    call HIST_reg( AD_HIST_id(I_TKE_RS   ), 'TKE_RS',    'resolved scale TKE',             'm2/s2',  ndim=3 )

    call HIST_reg( AD_HIST_id(I_ENGT)     , 'ENGT',      'total energy',                   'J/m3',   ndim=3 )
    call HIST_reg( AD_HIST_id(I_ENGP)     , 'ENGP',      'potential energy',               'J/m3',   ndim=3 )
    call HIST_reg( AD_HIST_id(I_ENGK)     , 'ENGK',      'kinetic energy',                 'J/m3',   ndim=3 )
    call HIST_reg( AD_HIST_id(I_ENGI)     , 'ENGI',      'internal energy',                'J/m3',   ndim=3 )

    !-----< monitor output setup >-----

    call MONIT_reg( AD_MONIT_id(I_QDRY)        , 'QDRY',         'dry air mass',      'kg', ndim=3, isflux=.false. )
    call MONIT_reg( AD_MONIT_id(I_QTOT)        , 'QTOT',         'water mass',        'kg', ndim=3, isflux=.false. )
    call MONIT_reg( AD_MONIT_id(I_EVAP)        , 'EVAP',         'evaporation',       'kg', ndim=2, isflux=.true.  )
    call MONIT_reg( AD_MONIT_id(I_PRCP)        , 'PRCP',         'precipitation',     'kg', ndim=2, isflux=.true.  )

    call MONIT_reg( AD_MONIT_id(I_ENGT)        , 'ENGT',         'total     energy',  'J',  ndim=3, isflux=.false. )
    call MONIT_reg( AD_MONIT_id(I_ENGP)        , 'ENGP',         'potential energy',  'J',  ndim=3, isflux=.false. )
    call MONIT_reg( AD_MONIT_id(I_ENGK)        , 'ENGK',         'kinetic   energy',  'J',  ndim=3, isflux=.false. )
    call MONIT_reg( AD_MONIT_id(I_ENGI)        , 'ENGI',         'internal  energy',  'J',  ndim=3, isflux=.false. )

    call MONIT_reg( AD_MONIT_id(I_ENGFLXT)     , 'ENGFLXT',      'total energy flux', 'J',  ndim=2, isflux=.true.  )

    call MONIT_reg( AD_MONIT_id(I_ENGSFC_SH)   , 'ENGSFC_SH',    'total     energy',  'J',  ndim=2, isflux=.true.  )
    call MONIT_reg( AD_MONIT_id(I_ENGSFC_LH)   , 'ENGSFC_LH',    'potential energy',  'J',  ndim=2, isflux=.true.  )
    call MONIT_reg( AD_MONIT_id(I_ENGSFC_RD)   , 'ENGSFC_RD',    'kinetic   energy',  'J',  ndim=2, isflux=.true.  )
    call MONIT_reg( AD_MONIT_id(I_ENGTOA_RD)   , 'ENGTOA_RD',    'internal  energy',  'J',  ndim=2, isflux=.true.  )

    call MONIT_reg( AD_MONIT_id(I_ENGSFC_LW_up), 'ENGSFC_LW_up', 'total     energy',  'J',  ndim=2, isflux=.true.  )
    call MONIT_reg( AD_MONIT_id(I_ENGSFC_LW_dn), 'ENGSFC_LW_dn', 'potential energy',  'J',  ndim=2, isflux=.true.  )
    call MONIT_reg( AD_MONIT_id(I_ENGSFC_SW_up), 'ENGSFC_SW_up', 'kinetic   energy',  'J',  ndim=2, isflux=.true.  )
    call MONIT_reg( AD_MONIT_id(I_ENGSFC_SW_dn), 'ENGSFC_SW_dn', 'internal  energy',  'J',  ndim=2, isflux=.true.  )

    call MONIT_reg( AD_MONIT_id(I_ENGTOA_LW_up), 'ENGTOA_LW_up', 'total     energy',  'J',  ndim=2, isflux=.true.  )
    call MONIT_reg( AD_MONIT_id(I_ENGTOA_LW_dn), 'ENGTOA_LW_dn', 'potential energy',  'J',  ndim=2, isflux=.true.  )
    call MONIT_reg( AD_MONIT_id(I_ENGTOA_SW_up), 'ENGTOA_SW_up', 'kinetic   energy',  'J',  ndim=2, isflux=.true.  )
    call MONIT_reg( AD_MONIT_id(I_ENGTOA_SW_dn), 'ENGTOA_SW_dn', 'internal  energy',  'J',  ndim=2, isflux=.true.  )

    if ( AD_HIST_id(I_W) > 0 ) then
       AD_PREP_sw(I_W) = 1
    endif
    if ( AD_HIST_id(I_U) > 0 ) then
       AD_PREP_sw(I_U) = 1
    endif
    if ( AD_HIST_id(I_V) > 0 ) then
       AD_PREP_sw(I_V) = 1
    endif
    if ( AD_HIST_id(I_POTT) > 0 ) then
       AD_PREP_sw(I_POTT) = 1
    endif

    if (      AD_HIST_id(I_QDRY) > 0 &
         .OR. AD_MONIT_id(I_QDRY) > 0 ) then
       AD_PREP_sw(I_QDRY) = 1
    endif
    if (      AD_HIST_id (I_QTOT) > 0 &
         .OR. AD_MONIT_id(I_QTOT) > 0 ) then
       AD_PREP_sw(I_QDRY) = 1
       AD_PREP_sw(I_QTOT) = 1
    endif
    if ( AD_HIST_id(I_QHYD) > 0 ) then
       AD_PREP_sw(I_QHYD) = 1
    endif
    if ( AD_HIST_id(I_QLIQ) > 0 ) then
       AD_PREP_sw(I_QLIQ) = 1
    endif
    if ( AD_HIST_id(I_QICE) > 0 ) then
       AD_PREP_sw(I_QICE) = 1
    endif

    if ( AD_HIST_id(I_LWP)  > 0 ) then
       AD_PREP_sw(I_QLIQ) = 1
       AD_PREP_sw(I_LWP)  = 1
    endif
    if ( AD_HIST_id(I_IWP)  > 0 ) then
       AD_PREP_sw(I_QICE) = 1
       AD_PREP_sw(I_IWP)  = 1
    endif
    if ( AD_HIST_id(I_PW)  > 0 ) then
       AD_PREP_sw(I_PW)  = 1
    endif

    if ( AD_HIST_id(I_RTOT) > 0 ) then
       AD_PREP_sw(I_QDRY) = 1
       AD_PREP_sw(I_RTOT) = 1
    endif
    if ( AD_HIST_id(I_CPTOT) > 0 ) then
       AD_PREP_sw(I_QDRY)  = 1
       AD_PREP_sw(I_CPTOT) = 1
    endif
    if ( AD_HIST_id(I_PRES) > 0 ) then
       AD_PREP_sw(I_QDRY)  = 1
       AD_PREP_sw(I_RTOT)  = 1
       AD_PREP_sw(I_CPTOT) = 1
       AD_PREP_sw(I_PRES)  = 1
    endif
    if ( AD_HIST_id(I_TEMP) > 0 ) then
       AD_PREP_sw(I_QDRY)  = 1
       AD_PREP_sw(I_RTOT)  = 1
       AD_PREP_sw(I_CPTOT) = 1
       AD_PREP_sw(I_PRES)  = 1
       AD_PREP_sw(I_TEMP)  = 1
    endif

    if ( AD_HIST_id(I_POTL) > 0 ) then
       AD_PREP_sw(I_POTT)  = 1
       AD_PREP_sw(I_QDRY)  = 1
       AD_PREP_sw(I_RTOT)  = 1
       AD_PREP_sw(I_CPTOT) = 1
       AD_PREP_sw(I_PRES)  = 1
       AD_PREP_sw(I_TEMP)  = 1
       AD_PREP_sw(I_POTL)  = 1
    endif
    if (      AD_HIST_id(I_RHA) > 0 &
         .OR. AD_HIST_id(I_RHL) > 0 &
         .OR. AD_HIST_id(I_RHI) > 0 ) then
       AD_PREP_sw(I_QDRY)  = 1
       AD_PREP_sw(I_RTOT)  = 1
       AD_PREP_sw(I_CPTOT) = 1
       AD_PREP_sw(I_PRES)  = 1
       AD_PREP_sw(I_TEMP)  = 1
       AD_PREP_sw(I_QSAT)  = 1
    endif


    if ( AD_HIST_id (I_VOR) > 0 ) then
       AD_PREP_sw(I_VOR) = 1
    endif

    if ( AD_PREP_sw(I_DIV) > 0 ) then
       AD_PREP_sw(I_HDIV) = 1
    endif

    if ( AD_HIST_id(I_Uabs) > 0 ) then
       AD_PREP_sw(I_U)    = 1
       AD_PREP_sw(I_V)    = 1
       AD_PREP_sw(I_Uabs) = 1
    endif

    if (      AD_HIST_id(I_CAPE) > 0 &
         .OR. AD_HIST_id(I_CIN)  > 0 &
         .OR. AD_HIST_id(I_LCL)  > 0 &
         .OR. AD_HIST_id(I_LFC)  > 0 &
         .OR. AD_HIST_id(I_LNB)  > 0 ) then
       AD_PREP_sw(I_CAPE) = 1
       AD_PREP_sw(I_CIN)  = 1
       AD_PREP_sw(I_LCL)  = 1
       AD_PREP_sw(I_LFC)  = 1
       AD_PREP_sw(I_LNB)  = 1
    endif

    if ( AD_HIST_id(I_PBLH) > 0 ) then
       AD_PREP_sw(I_POTT) = 1
       AD_PREP_sw(I_PBLH) = 1
    endif

    if ( AD_HIST_id(I_MSE) > 0 ) then
       AD_PREP_sw(I_CPTOT) = 1
       AD_PREP_sw(I_TEMP)  = 1
       AD_PREP_sw(I_MSE)   = 1
    endif

    if ( AD_HIST_id(I_DENS_PRIM) > 0 ) then
       AD_PREP_sw(I_DENS_PRIM) = 1
       AD_PREP_sw(I_DENS_MEAN) = 1
    endif

    if ( AD_HIST_id(I_W_PRIM) > 0 ) then
       AD_PREP_sw(I_W)      = 1
       AD_PREP_sw(I_W_PRIM) = 1
       AD_PREP_sw(I_DENS_MEAN) = 1
       AD_PREP_sw(I_W_MEAN) = 1
    endif

    if ( AD_HIST_id(I_U_PRIM) > 0 ) then
       AD_PREP_sw(I_U)      = 1
       AD_PREP_sw(I_U_PRIM) = 1
       AD_PREP_sw(I_DENS_MEAN) = 1
       AD_PREP_sw(I_U_MEAN) = 1
    endif

    if ( AD_HIST_id(I_V_PRIM) > 0 ) then
       AD_PREP_sw(I_V)      = 1
       AD_PREP_sw(I_V_PRIM) = 1
       AD_PREP_sw(I_DENS_MEAN) = 1
       AD_PREP_sw(I_V_MEAN) = 1
    endif

    if ( AD_HIST_id(I_POTT_PRIM) > 0 ) then
       AD_PREP_sw(I_POTT)      = 1
       AD_PREP_sw(I_POTT_PRIM) = 1
       AD_PREP_sw(I_DENS_MEAN) = 1
       AD_PREP_sw(I_POTT_MEAN) = 1
    endif

    if ( AD_HIST_id(I_DENS_MEAN) > 0 ) then
       AD_PREP_sw(I_DENS_MEAN) = 1
    endif

    if ( AD_HIST_id(I_W_MEAN) > 0 ) then
       AD_PREP_sw(I_W_MEAN) = 1
    endif

    if ( AD_HIST_id(I_U_MEAN) > 0 ) then
       AD_PREP_sw(I_U_MEAN) = 1
    endif

    if ( AD_HIST_id(I_V_MEAN) > 0 ) then
       AD_PREP_sw(I_V_MEAN) = 1
    endif

    if ( AD_HIST_id(I_POTT_MEAN) > 0 ) then
       AD_PREP_sw(I_POTT_MEAN) = 1
    end if

    if ( AD_HIST_id(I_T_MEAN) > 0 ) then
       AD_PREP_sw(I_T_MEAN) = 1
    end if

    if ( AD_HIST_id(I_QV_MEAN) > 0 ) then
       AD_PREP_sw(I_QV_MEAN) = 1
    end if

    if ( AD_HIST_id(I_QHYD_MEAN) > 0 ) then
       AD_PREP_sw(I_QHYD) = 1
       AD_PREP_sw(I_QHYD_MEAN) = 1
    end if

    if ( AD_HIST_id(I_QLIQ_MEAN) > 0 ) then
       AD_PREP_sw(I_QLIQ) = 1
       AD_PREP_sw(I_QLIQ_MEAN) = 1
    end if

    if ( AD_HIST_id(I_QICE_MEAN) > 0 ) then
       AD_PREP_sw(I_QICE) = 1
       AD_PREP_sw(I_QICE_MEAN) = 1
    end if

    if ( AD_HIST_id(I_W_PRIM2) > 0 ) then
       AD_PREP_sw(I_W)       = 1
       AD_PREP_sw(I_W_PRIM)  = 1
       AD_PREP_sw(I_W_PRIM2) = 1
       AD_PREP_sw(I_DENS_MEAN) = 1
       AD_PREP_sw(I_W_MEAN)    = 1
    endif

    if ( AD_HIST_id(I_PT_W_PRIM) > 0 ) then
       AD_PREP_sw(I_W)         = 1
       AD_PREP_sw(I_W_PRIM)    = 1
       AD_PREP_sw(I_POTT)      = 1
       AD_PREP_sw(I_POTT_PRIM) = 1
       AD_PREP_sw(I_PT_W_PRIM) = 1
       AD_PREP_sw(I_DENS_MEAN) = 1
       AD_PREP_sw(I_W_MEAN)    = 1
       AD_PREP_sw(I_POTT_MEAN) = 1
    endif

    if ( AD_HIST_id(I_W_PRIM3) > 0 ) then
       AD_PREP_sw(I_W)       = 1
       AD_PREP_sw(I_W_PRIM)  = 1
       AD_PREP_sw(I_W_PRIM3) = 1
       AD_PREP_sw(I_DENS_MEAN) = 1
       AD_PREP_sw(I_W_MEAN)  = 1
    endif

    if ( AD_HIST_id(I_TKE_RS) > 0 ) then
       AD_PREP_sw(I_W)      = 1
       AD_PREP_sw(I_U)      = 1
       AD_PREP_sw(I_V)      = 1
       AD_PREP_sw(I_W_PRIM) = 1
       AD_PREP_sw(I_U_PRIM) = 1
       AD_PREP_sw(I_V_PRIM) = 1
       AD_PREP_sw(I_TKE_RS) = 1
       AD_PREP_sw(I_DENS_MEAN) = 1
       AD_PREP_sw(I_W_MEAN) = 1
       AD_PREP_sw(I_U_MEAN) = 1
       AD_PREP_sw(I_V_MEAN) = 1
    endif

    if (      AD_HIST_id (I_ENGP) > 0 &
         .OR. AD_MONIT_id(I_ENGP) > 0 ) then
       AD_PREP_sw(I_ENGP) = 1
    endif
    if (      AD_HIST_id (I_ENGK) > 0 &
         .OR. AD_MONIT_id(I_ENGK) > 0 ) then
       AD_PREP_sw(I_W) = 1
       AD_PREP_sw(I_U) = 1
       AD_PREP_sw(I_V) = 1
       AD_PREP_sw(I_ENGK) = 1
    endif
    if (      AD_HIST_id (I_ENGI) > 0 &
         .OR. AD_MONIT_id(I_ENGI) > 0 ) then
       AD_PREP_sw(I_QDRY)  = 1
       AD_PREP_sw(I_RTOT)  = 1
       AD_PREP_sw(I_CPTOT) = 1
       AD_PREP_sw(I_PRES)  = 1
       AD_PREP_sw(I_TEMP)  = 1
       AD_PREP_sw(I_ENGI)  = 1
    endif
    if (      AD_HIST_id (I_ENGT) > 0 &
         .OR. AD_MONIT_id(I_ENGT) > 0 ) then
       AD_PREP_sw(I_ENGP)  = 1
       AD_PREP_sw(I_W)  = 1
       AD_PREP_sw(I_U)  = 1
       AD_PREP_sw(I_V)  = 1
       AD_PREP_sw(I_ENGK)  = 1
       AD_PREP_sw(I_QDRY)  = 1
       AD_PREP_sw(I_RTOT)  = 1
       AD_PREP_sw(I_CPTOT) = 1
       AD_PREP_sw(I_PRES)  = 1
       AD_PREP_sw(I_TEMP)  = 1
       AD_PREP_sw(I_ENGI)  = 1
       AD_PREP_sw(I_ENGT)  = 1
    endif

    return
  end subroutine ATMOS_vars_setup

  !-----------------------------------------------------------------------------
  !> HALO Communication
  subroutine ATMOS_vars_fillhalo( &
       FILL_BND )
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    logical, intent(in), optional :: FILL_BND

    logical :: FILL_BND_
    integer :: i, j, iq
    !---------------------------------------------------------------------------

    FILL_BND_ = .false.
    if ( present(FILL_BND) ) FILL_BND_ = FILL_BND

    !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
    do j  = JSB, JEB
    do i  = ISB, IEB
       DENS(   1:KS-1,i,j) = DENS(KS,i,j)
       MOMZ(   1:KS-1,i,j) = MOMZ(KS,i,j)
       MOMX(   1:KS-1,i,j) = MOMX(KS,i,j)
       MOMY(   1:KS-1,i,j) = MOMY(KS,i,j)
       RHOT(   1:KS-1,i,j) = RHOT(KS,i,j)
       DENS(KE+1:KA,  i,j) = DENS(KE,i,j)
       MOMZ(KE+1:KA,  i,j) = MOMZ(KE,i,j)
       MOMX(KE+1:KA,  i,j) = MOMX(KE,i,j)
       MOMY(KE+1:KA,  i,j) = MOMY(KE,i,j)
       RHOT(KE+1:KA,  i,j) = RHOT(KE,i,j)
    enddo
    enddo

    !$omp parallel do private(i,j,iq) OMP_SCHEDULE_ collapse(3)
    do iq = 1, QA
    do j  = JSB, JEB
    do i  = ISB, IEB
       QTRC(   1:KS-1,i,j,iq) = QTRC(KS,i,j,iq)
       QTRC(KE+1:KA,  i,j,iq) = QTRC(KE,i,j,iq)
    enddo
    enddo
    enddo

    call COMM_vars8( DENS(:,:,:), 1 )
    call COMM_vars8( MOMZ(:,:,:), 2 )
    call COMM_vars8( MOMX(:,:,:), 3 )
    call COMM_vars8( MOMY(:,:,:), 4 )
    call COMM_vars8( RHOT(:,:,:), 5 )
    call COMM_wait ( DENS(:,:,:), 1, FILL_BND_ )
    call COMM_wait ( MOMZ(:,:,:), 2, FILL_BND_ )
    call COMM_wait ( MOMX(:,:,:), 3, FILL_BND_ )
    call COMM_wait ( MOMY(:,:,:), 4, FILL_BND_ )
    call COMM_wait ( RHOT(:,:,:), 5, FILL_BND_ )

    do iq = 1, QA
       call COMM_vars8( QTRC(:,:,:,iq), iq )
    enddo
    do iq = 1, QA
       call COMM_wait ( QTRC(:,:,:,iq), iq, FILL_BND_ )
    enddo

    return
  end subroutine ATMOS_vars_fillhalo

  !-----------------------------------------------------------------------------
  !> Open restart file for reading atmospheric variables
  subroutine ATMOS_vars_restart_open
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       GRAV  => CONST_GRAV
    use scale_time, only: &
       TIME_gettimelabel
    use scale_fileio, only: &
       FILEIO_open
    use scale_time, only: &
       TIME_gettimelabel
    use scale_atmos_thermodyn, only: &
       THERMODYN_qd        => ATMOS_THERMODYN_qd,        &
       THERMODYN_temp_pres => ATMOS_THERMODYN_temp_pres
    use mod_atmos_admin, only: &
       ATMOS_USE_AVERAGE, &
       ATMOS_sw_dyn,      &
       ATMOS_sw_phy_mp,   &
       ATMOS_sw_phy_ae,   &
       ATMOS_sw_phy_ch,   &
       ATMOS_sw_phy_rd,   &
       ATMOS_sw_phy_sf,   &
       ATMOS_sw_phy_tb,   &
       ATMOS_sw_phy_cp
    use mod_atmos_dyn_vars, only: &
       ATMOS_DYN_vars_restart_open
    use mod_atmos_phy_mp_vars, only: &
       ATMOS_PHY_MP_vars_restart_open
    use mod_atmos_phy_ae_vars, only: &
       ATMOS_PHY_AE_vars_restart_open
    use mod_atmos_phy_ch_vars, only: &
       ATMOS_PHY_CH_vars_restart_open
    use mod_atmos_phy_rd_vars, only: &
       ATMOS_PHY_RD_vars_restart_open
    use mod_atmos_phy_sf_vars, only: &
       ATMOS_PHY_SF_vars_restart_open
    use mod_atmos_phy_tb_vars, only: &
       ATMOS_PHY_TB_vars_restart_open
    use mod_atmos_phy_cp_vars, only: &
       ATMOS_PHY_CP_vars_restart_open
    implicit none

    character(len=20)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input restart file (ATMOS) ***'

    if ( ATMOS_RESTART_IN_BASENAME /= '' ) then

       if ( ATMOS_RESTART_IN_POSTFIX_TIMELABEL ) then
          call TIME_gettimelabel( timelabel )
          basename = trim(ATMOS_RESTART_IN_BASENAME)//'_'//trim(timelabel)
       else
          basename = trim(ATMOS_RESTART_IN_BASENAME)
       endif

       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(basename)

       call FILEIO_open( restart_fid, basename )
    else
       write(*,*) '*** restart file for atmosphere is not specified. STOP!'
       call PRC_MPIstop
    endif

    if ( ATMOS_USE_AVERAGE ) then
       DENS_av(:,:,:)   = DENS(:,:,:)
       MOMZ_av(:,:,:)   = MOMZ(:,:,:)
       MOMX_av(:,:,:)   = MOMX(:,:,:)
       MOMY_av(:,:,:)   = MOMY(:,:,:)
       RHOT_av(:,:,:)   = RHOT(:,:,:)
       QTRC_av(:,:,:,:) = QTRC(:,:,:,:)
    endif

    if( ATMOS_sw_dyn )    call ATMOS_DYN_vars_restart_open
    if( ATMOS_sw_phy_mp ) call ATMOS_PHY_MP_vars_restart_open
    if( ATMOS_sw_phy_ae ) call ATMOS_PHY_AE_vars_restart_open
    if( ATMOS_sw_phy_ch ) call ATMOS_PHY_CH_vars_restart_open
    if( ATMOS_sw_phy_rd ) call ATMOS_PHY_RD_vars_restart_open
    if( ATMOS_sw_phy_sf ) call ATMOS_PHY_SF_vars_restart_open
    if( ATMOS_sw_phy_tb ) call ATMOS_PHY_TB_vars_restart_open
    if( ATMOS_sw_phy_cp ) call ATMOS_PHY_CP_vars_restart_open

    return
  end subroutine ATMOS_vars_restart_open

  !-----------------------------------------------------------------------------
  !> Read restart of atmospheric variables
  subroutine ATMOS_vars_restart_read
    use scale_process, only: &
       PRC_MPIstop
    use scale_time, only: &
       TIME_gettimelabel
    use scale_fileio, only: &
       FILEIO_read, &
       FILEIO_flush
    use scale_atmos_thermodyn, only: &
       THERMODYN_qd        => ATMOS_THERMODYN_qd,        &
       THERMODYN_temp_pres => ATMOS_THERMODYN_temp_pres
    use mod_atmos_admin, only: &
       ATMOS_USE_AVERAGE, &
       ATMOS_sw_dyn,      &
       ATMOS_sw_phy_mp,   &
       ATMOS_sw_phy_ae,   &
       ATMOS_sw_phy_ch,   &
       ATMOS_sw_phy_rd,   &
       ATMOS_sw_phy_sf,   &
       ATMOS_sw_phy_tb,   &
       ATMOS_sw_phy_cp
    use mod_atmos_dyn_vars, only: &
       ATMOS_DYN_vars_restart_read
    use mod_atmos_phy_mp_vars, only: &
       ATMOS_PHY_MP_vars_restart_read
    use mod_atmos_phy_ae_vars, only: &
       ATMOS_PHY_AE_vars_restart_read
    use mod_atmos_phy_ch_vars, only: &
       ATMOS_PHY_CH_vars_restart_read
    use mod_atmos_phy_rd_vars, only: &
       ATMOS_PHY_RD_vars_restart_read
    use mod_atmos_phy_sf_vars, only: &
       ATMOS_PHY_SF_vars_restart_read
    use mod_atmos_phy_tb_vars, only: &
       ATMOS_PHY_TB_vars_restart_read
    use mod_atmos_phy_cp_vars, only: &
       ATMOS_PHY_CP_vars_restart_read
    implicit none

    integer  :: i, j, iq
    !---------------------------------------------------------------------------

    if ( restart_fid .NE. -1 ) then

       call FILEIO_read( DENS(:,:,:),                            & ! [OUT]
                         restart_fid, VAR_NAME(1), 'ZXY', step=1 ) ! [IN]
       call FILEIO_read( MOMZ(:,:,:),                            & ! [OUT]
                         restart_fid, VAR_NAME(2), 'ZXY', step=1 ) ! [IN]
       call FILEIO_read( MOMX(:,:,:),                            & ! [OUT]
                         restart_fid, VAR_NAME(3), 'ZXY', step=1 ) ! [IN]
       call FILEIO_read( MOMY(:,:,:),                            & ! [OUT]
                         restart_fid, VAR_NAME(4), 'ZXY', step=1 ) ! [IN]
       call FILEIO_read( RHOT(:,:,:),                            & ! [OUT]
                         restart_fid, VAR_NAME(5), 'ZXY', step=1 ) ! [IN]

       do iq = 1, QA
          call FILEIO_read( QTRC(:,:,:,iq),                             & ! [OUT]
                            restart_fid, TRACER_NAME(iq), 'ZXY', step=1 ) ! [IN]
       enddo

       if ( IO_AGGREGATE ) then
          call FILEIO_flush( restart_fid )
          ! X/Y halos have been read from file

          ! fill k halos
          do j  = 1, JA
          do i  = 1, IA
             DENS(   1:KS-1,i,j) = DENS(KS,i,j)
             MOMZ(   1:KS-1,i,j) = MOMZ(KS,i,j)
             MOMX(   1:KS-1,i,j) = MOMX(KS,i,j)
             MOMY(   1:KS-1,i,j) = MOMY(KS,i,j)
             RHOT(   1:KS-1,i,j) = RHOT(KS,i,j)
             DENS(KE+1:KA,  i,j) = DENS(KE,i,j)
             MOMZ(KE+1:KA,  i,j) = MOMZ(KE,i,j)
             MOMX(KE+1:KA,  i,j) = MOMX(KE,i,j)
             MOMY(KE+1:KA,  i,j) = MOMY(KE,i,j)
             RHOT(KE+1:KA,  i,j) = RHOT(KE,i,j)
          enddo
          enddo
       else
          call ATMOS_vars_fillhalo
       end if

       call ATMOS_vars_total
    else
       write(*,*) '*** invalid restart file ID for atmosphere. STOP!'
       call PRC_MPIstop
    endif

    if ( ATMOS_USE_AVERAGE ) then
       DENS_av(:,:,:)   = DENS(:,:,:)
       MOMZ_av(:,:,:)   = MOMZ(:,:,:)
       MOMX_av(:,:,:)   = MOMX(:,:,:)
       MOMY_av(:,:,:)   = MOMY(:,:,:)
       RHOT_av(:,:,:)   = RHOT(:,:,:)
       QTRC_av(:,:,:,:) = QTRC(:,:,:,:)
    endif

    if ( ATMOS_sw_dyn )    call ATMOS_DYN_vars_restart_read
    if ( ATMOS_sw_phy_mp ) call ATMOS_PHY_MP_vars_restart_read
    if ( ATMOS_sw_phy_ae ) call ATMOS_PHY_AE_vars_restart_read
    if ( ATMOS_sw_phy_ch ) call ATMOS_PHY_CH_vars_restart_read
    if ( ATMOS_sw_phy_rd ) call ATMOS_PHY_RD_vars_restart_read
    if ( ATMOS_sw_phy_sf ) call ATMOS_PHY_SF_vars_restart_read
    if ( ATMOS_sw_phy_tb ) call ATMOS_PHY_TB_vars_restart_read
    if ( ATMOS_sw_phy_cp ) call ATMOS_PHY_CP_vars_restart_read

    return
  end subroutine ATMOS_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Set pressure for history output
  subroutine ATMOS_vars_history_setpres
    use scale_grid_real, only: &
       REAL_CZ, &
       REAL_Z1
    use scale_topography, only: &
       TOPO_Zsfc
    use scale_atmos_bottom, only: &
       BOTTOM_estimate => ATMOS_BOTTOM_estimate
    use scale_history, only: &
       HIST_setpres
    implicit none

    real(RP) :: SFC_DENS(IA,JA)
    real(RP) :: SFC_PRES(IA,JA)
    !---------------------------------------------------------------------------

    call BOTTOM_estimate( DENS     (:,:,:), & ! [IN]
                          PRES     (:,:,:), & ! [IN]
                          REAL_CZ  (:,:,:), & ! [IN]
                          TOPO_Zsfc(:,:),   & ! [IN]
                          REAL_Z1  (:,:),   & ! [IN]
                          SFC_DENS (:,:),   & ! [OUT]
                          SFC_PRES (:,:)    ) ! [OUT]

    call HIST_setpres( PRES    (:,:,:),  & ! [IN]
                       SFC_PRES(:,:)     ) ! [IN]

    return
  end subroutine ATMOS_vars_history_setpres

  !-----------------------------------------------------------------------------
  !> Check and compare between last data and sample data
  subroutine ATMOS_vars_restart_check
    use scale_time, only: &
       TIME_gettimelabel
    use scale_process, only: &
       PRC_myrank
    use scale_fileio, only: &
       FILEIO_open, &
       FILEIO_read, &
       FILEIO_flush, &
       FILEIO_close
    implicit none

    real(RP) :: DENS_check(KA,IA,JA)    ! Density    [kg/m3]
    real(RP) :: MOMZ_check(KA,IA,JA)    ! momentum z [kg/s/m2]
    real(RP) :: MOMX_check(KA,IA,JA)    ! momentum x [kg/s/m2]
    real(RP) :: MOMY_check(KA,IA,JA)    ! momentum y [kg/s/m2]
    real(RP) :: RHOT_check(KA,IA,JA)    ! DENS * POTT [K*kg/m3]
    real(RP) :: QTRC_check(KA,IA,JA,QA) ! tracer mixing ratio [kg/kg]

    character(len=H_LONG) :: basename
    character(len=20)     :: timelabel

    logical :: datacheck
    integer :: k, i, j, iq
    integer :: fid
    !---------------------------------------------------------------------------

    call PROF_rapstart('Debug')

    write(*,*) 'Compare last Data with ', trim(ATMOS_RESTART_CHECK_BASENAME), 'on PE=', PRC_myrank
    write(*,*) '*** criterion = ', ATMOS_RESTART_CHECK_CRITERION
    datacheck = .true.

    basename = ATMOS_RESTART_CHECK_BASENAME

    call FILEIO_open( fid, basename )

    call FILEIO_read( DENS_check(:,:,:), fid, 'DENS', 'ZXY', step=1 )
    call FILEIO_read( MOMZ_check(:,:,:), fid, 'MOMZ', 'ZXY', step=1 )
    call FILEIO_read( MOMX_check(:,:,:), fid, 'MOMX', 'ZXY', step=1 )
    call FILEIO_read( MOMY_check(:,:,:), fid, 'MOMY', 'ZXY', step=1 )
    call FILEIO_read( RHOT_check(:,:,:), fid, 'RHOT', 'ZXY', step=1 )
    do iq = 1, QA
       call FILEIO_read( QTRC_check(:,:,:,iq), fid, TRACER_NAME(iq), 'ZXY', step=1 )
    end do
    if ( IO_AGGREGATE ) call FILEIO_flush( fid )

    call FILEIO_close( fid ) ! [IN]

    do k = KS, KE
    do j = JS, JE
    do i = IS, IE
       if ( abs( DENS(k,i,j)-DENS_check(k,i,j) ) > ATMOS_RESTART_CHECK_CRITERION ) then
          write(*,*) 'xxx there is the difference  : ', DENS(k,i,j)-DENS_check(k,i,j)
          write(*,*) 'xxx at (PE-id,k,i,j,varname) : ', PRC_myrank, k, i, j, 'DENS'
          datacheck = .false.
       endif
    enddo
    enddo
    enddo

    do k = KS, KE
    do j = JS, JE
    do i = IS, IE
       if ( abs( MOMZ(k,i,j)-MOMZ_check(k,i,j) ) > ATMOS_RESTART_CHECK_CRITERION ) then
          write(*,*) 'xxx there is the difference  : ', MOMZ(k,i,j)-MOMZ_check(k,i,j)
          write(*,*) 'xxx at (PE-id,k,i,j,varname) : ', PRC_myrank, k, i, j, 'MOMZ'
          datacheck = .false.
       endif
    enddo
    enddo
    enddo

    do k = KS, KE
    do j = JS, JE
    do i = IS, IE
       if ( abs( MOMX(k,i,j)-MOMX_check(k,i,j) ) > ATMOS_RESTART_CHECK_CRITERION ) then
          write(*,*) 'xxx there is the difference  : ', MOMX(k,i,j)-MOMX_check(k,i,j)
          write(*,*) 'xxx at (PE-id,k,i,j,varname) : ', PRC_myrank, k, i, j, 'MOMX'
          datacheck = .false.
       endif
    enddo
    enddo
    enddo

    do k = KS, KE
    do j = JS, JE
    do i = IS, IE
       if ( abs( MOMY(k,i,j)-MOMY_check(k,i,j) ) > ATMOS_RESTART_CHECK_CRITERION ) then
          write(*,*) 'xxx there is the difference  : ', MOMY(k,i,j)-MOMY_check(k,i,j)
          write(*,*) 'xxx at (PE-id,k,i,j,varname) : ', PRC_myrank, k, i, j, 'MOMY'
          datacheck = .false.
       endif
    enddo
    enddo
    enddo

    do k = KS, KE
    do j = JS, JE
    do i = IS, IE
       if ( abs( RHOT(k,i,j)-RHOT_check(k,i,j) ) > ATMOS_RESTART_CHECK_CRITERION ) then
          write(*,*) 'xxx there is the difference  : ', RHOT(k,i,j)-RHOT_check(k,i,j)
          write(*,*) 'xxx at (PE-id,k,i,j,varname) : ', PRC_myrank, k, i, j, 'RHOT'
          datacheck = .false.
       endif
    enddo
    enddo
    enddo

    do iq = 1, QA
       do k = KS, KE
       do j = JS, JE
       do i = IS, IE
          if ( abs( QTRC(k,i,j,iq)-QTRC_check(k,i,j,iq) ) > ATMOS_RESTART_CHECK_CRITERION ) then
             write(*,*) 'xxx there is the difference  : ', QTRC(k,i,j,iq)-QTRC_check(k,i,j,iq)
             write(*,*) 'xxx at (PE-id,k,i,j,varname) : ', PRC_myrank, k, i, j, TRACER_NAME(iq)
             datacheck = .false.
          endif
       enddo
       enddo
       enddo
    enddo

    if (datacheck) then
       if( IO_L ) write(IO_FID_LOG,*) 'Data Check Clear.'
       write(*,*) 'Data Check Clear.'
    else
       if( IO_L ) write(IO_FID_LOG,*) 'Data Check Failed. See std. output.'
       write(*,*) 'Data Check Failed.'
    endif

    call PROF_rapend('Debug')

    return
  end subroutine ATMOS_vars_restart_check

  !-----------------------------------------------------------------------------
  !> History output set for atmospheric variables
  subroutine ATMOS_vars_history
    use scale_const, only: &
       GRAV  => CONST_GRAV,  &
       Rdry  => CONST_Rdry,  &
       Rvap  => CONST_Rvap,  &
       CPdry => CONST_CPdry, &
       CVdry => CONST_CVdry, &
       P00   => CONST_PRE00
    use scale_grid, only: &
       RCDX => GRID_RCDX, &
       RCDY => GRID_RCDY
    use scale_grid_real, only: &
       REAL_CZ, &
       REAL_FZ
    use scale_comm, only: &
       COMM_horizontal_mean
    use scale_history, only: &
       HIST_in
    use scale_atmos_thermodyn, only: &
       THERMODYN_qd => ATMOS_THERMODYN_qd
    use scale_atmos_hydrometeor, only: &
       HYDROMETEOR_LHV => ATMOS_HYDROMETEOR_LHV, &
       I_QV, &
       I_QC, &
       QHS,  &
       QHE,  &
       QLS,  &
       QLE,  &
       QIS,  &
       QIE,  &
       LHV,  &
       LHF
    use scale_atmos_saturation, only: &
       SATURATION_psat_all => ATMOS_SATURATION_psat_all, &
       SATURATION_psat_liq => ATMOS_SATURATION_psat_liq, &
       SATURATION_psat_ice => ATMOS_SATURATION_psat_ice
    use scale_atmos_adiabat, only: &
       ADIABAT_cape => ATMOS_ADIABAT_cape
    use mod_atmos_phy_cp_vars, only: &
       SFLX_rain_CP => ATMOS_PHY_CP_SFLX_rain
    use mod_atmos_phy_mp_vars, only: &
       SFLX_rain_MP => ATMOS_PHY_MP_SFLX_rain, &
       SFLX_snow_MP => ATMOS_PHY_MP_SFLX_snow
    implicit none

    real(RP) :: QDRY  (KA,IA,JA) ! dry air            [kg/kg]
    real(RP) :: QTOT  (KA,IA,JA) ! total water        [kg/kg]
    real(RP) :: QHYD  (KA,IA,JA) ! total hydrometeor  [kg/kg]
    real(RP) :: QLIQ  (KA,IA,JA) ! total liquid water [kg/kg]
    real(RP) :: QICE  (KA,IA,JA) ! total ice water    [kg/kg]
    real(RP) :: RHOQ  (KA,IA,JA)

    real(RP) :: LWP   (IA,JA)    ! liquid water path  [g/m2]
    real(RP) :: IWP   (IA,JA)    ! ice    water path  [g/m2]
    real(RP) :: PW    (IA,JA)    ! precipitable water [g/m2]

    real(RP) :: RTOT  (KA,IA,JA) ! Total gas constant  [J/kg/K]
    real(RP) :: CPTOT (KA,IA,JA) ! Total heat capacity [J/kg/K]
    real(RP) :: CVTOT (KA,IA,JA) ! Total heat capacity [J/kg/K]
    real(RP) :: CPovCV(KA,IA,JA) ! Cp/Cv

    real(RP) :: POTL  (KA,IA,JA) ! liquid water potential temperature [K]
    real(RP) :: RHA   (KA,IA,JA) ! relative humidity (liquid+ice)      [%]
    real(RP) :: RHL   (KA,IA,JA) ! relative humidity against to liquid [%]
    real(RP) :: RHI   (KA,IA,JA) ! relative humidity against to ice    [%]

    real(RP) :: VOR   (KA,IA,JA) ! vertical vorticity    [1/s]
    real(RP) :: DIV   (KA,IA,JA) ! divergence            [1/s]
    real(RP) :: HDIV  (KA,IA,JA) ! horizontal divergence [1/s]
    real(RP) :: Uabs  (KA,IA,JA) ! absolute velocity     [m/s]

    real(RP) :: CAPE  (IA,JA)    ! CAPE [m2/s2]
    real(RP) :: CIN   (IA,JA)    ! CIN [m2/s2]
    real(RP) :: LCL   (IA,JA)    ! LCL height [m]
    real(RP) :: LFC   (IA,JA)    ! LFC height [m]
    real(RP) :: LNB   (IA,JA)    ! LNB height [m]

    real(RP) :: PBLH  (IA,JA)    ! PBL height [m]
    real(RP) :: POTTv (KA,IA,JA) ! vertual potential temperature [K]
    real(RP) :: fact

    real(RP) :: MSE      (KA,IA,JA) ! MSE        [m2/s2]
    real(RP) :: LHV_local(KA,IA,JA) ! latent heat for vaporization [m2/s2]

    real(RP) :: PREC  (IA,JA)    ! surface precipitation rate CP+MP [kg/m2/s]

    real(RP) :: DENS_PRIM(KA,IA,JA) ! horiz. deviation of density    [kg/m3]
    real(RP) :: W_PRIM   (KA,IA,JA) ! horiz. deviation of w          [m/s]
    real(RP) :: U_PRIM   (KA,IA,JA) ! horiz. deviation of u          [m/s]
    real(RP) :: V_PRIM   (KA,IA,JA) ! horiz. deviation of v          [m/s]
    real(RP) :: POTT_PRIM(KA,IA,JA) ! horiz. deviation of pot. temp. [K]
    real(RP) :: W_PRIM2  (KA,IA,JA) ! variance of w                  [m2/s2]
    real(RP) :: PT_W_PRIM(KA,IA,JA) ! resolved scale heat flux       [W/s]
    real(RP) :: W_PRIM3  (KA,IA,JA) ! skewness of w                  [m3/s3]
    real(RP) :: TKE_RS   (KA,IA,JA) ! resolved scale TKE             [m2/s2]
    real(RP) :: DENS_MEAN(KA)       ! horiz. mean of density         [kg/m3]
    real(RP) :: W_MEAN   (KA)       ! horiz. mean of w               [m/s]
    real(RP) :: U_MEAN   (KA)       ! horiz. mean of u               [m/s]
    real(RP) :: V_MEAN   (KA)       ! horiz. mean of v               [m/s]
    real(RP) :: PT_MEAN  (KA)       ! horiz. mean of pot.            [K]
    real(RP) :: T_MEAN   (KA)       ! horiz. mean of t               [K]
    real(RP) :: QV_MEAN  (KA)       ! horiz. mean of QV
    real(RP) :: QHYD_MEAN(KA)       ! horiz. mean of QHYD
    real(RP) :: QLIQ_MEAN(KA)       ! horiz. mean of QLIQ
    real(RP) :: QICE_MEAN(KA)       ! horiz. mean of QICE

    real(RP) :: ENGT  (KA,IA,JA) ! total     energy [J/m3]
    real(RP) :: ENGP  (KA,IA,JA) ! potential energy [J/m3]
    real(RP) :: ENGK  (KA,IA,JA) ! kinetic   energy [J/m3]
    real(RP) :: ENGI  (KA,IA,JA) ! internal  energy [J/m3]

    real(RP) :: PSAT  (KA,IA,JA)
    real(RP) :: UH    (KA,IA,JA)
    real(RP) :: VH    (KA,IA,JA)

    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    ! value check for prognostic variables
    if ( ATMOS_VARS_CHECKRANGE ) then
       call VALCHECK( DENS(:,:,:),    0.0_RP,    2.0_RP, VAR_NAME(I_DENS), __FILE__, __LINE__ )
       call VALCHECK( MOMZ(:,:,:), -200.0_RP,  200.0_RP, VAR_NAME(I_MOMZ), __FILE__, __LINE__ )
       call VALCHECK( MOMX(:,:,:), -200.0_RP,  200.0_RP, VAR_NAME(I_MOMX), __FILE__, __LINE__ )
       call VALCHECK( MOMY(:,:,:), -200.0_RP,  200.0_RP, VAR_NAME(I_MOMY), __FILE__, __LINE__ )
       call VALCHECK( RHOT(:,:,:),    0.0_RP, 1000.0_RP, VAR_NAME(I_RHOT), __FILE__, __LINE__ )
    endif

    ! history output of prognostic variables
    call HIST_in( DENS(:,:,:), VAR_NAME(I_DENS), VAR_DESC(I_DENS), VAR_UNIT(I_DENS) )
    call HIST_in( MOMZ(:,:,:), VAR_NAME(I_MOMZ), VAR_DESC(I_MOMZ), VAR_UNIT(I_MOMZ) )
    call HIST_in( MOMX(:,:,:), VAR_NAME(I_MOMX), VAR_DESC(I_MOMX), VAR_UNIT(I_MOMX) )
    call HIST_in( MOMY(:,:,:), VAR_NAME(I_MOMY), VAR_DESC(I_MOMY), VAR_UNIT(I_MOMY) )
    call HIST_in( RHOT(:,:,:), VAR_NAME(I_RHOT), VAR_DESC(I_RHOT), VAR_UNIT(I_RHOT) )
    do iq = 1, QA
       call HIST_in( QTRC(:,:,:,iq), TRACER_NAME(iq), TRACER_DESC(iq), TRACER_UNIT(iq) )
    enddo

    ! prepare and history output of diagnostic variables

!    if ( AD_PREP_sw(I_W) > 0 ) then
!       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!       do j = JS, JE
!       do i = IS, IE
!       do k = KS, KE
!          W(k,i,j) = 0.5_RP * ( MOMZ(k-1,i,j)+MOMZ(k,i,j) ) / DENS(k,i,j)
!       enddo
!       enddo
!       enddo
!    endif
!
!    if ( AD_PREP_sw(I_U) > 0 ) then
!       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!       do j = JS, JE
!       do i = IS, IE
!       do k = KS, KE
!          U(k,i,j) = 0.5_RP * ( MOMX(k,i-1,j)+MOMX(k,i,j) ) / DENS(k,i,j)
!       enddo
!       enddo
!       enddo
!    endif
!
!    if ( AD_PREP_sw(I_V) > 0 ) then
!       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!       do j = JS, JE
!       do i = IS, IE
!       do k = KS, KE
!          V(k,i,j) = 0.5_RP * ( MOMY(k,i,j-1)+MOMY(k,i,j) ) / DENS(k,i,j)
!       enddo
!       enddo
!       enddo
!    endif
!
!    if ( AD_PREP_sw(I_POTT) > 0 ) then
!       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!       do j = JS, JE
!       do i = IS, IE
!       do k = KS, KE
!          POTT(k,i,j) = RHOT(k,i,j) / DENS(k,i,j)
!       enddo
!       enddo
!       enddo
!    endif

    if ( AD_PREP_sw(I_QDRY) > 0 ) then
       call THERMODYN_qd( QDRY(:,:,:),   & ! [OUT]
                          QTRC(:,:,:,:), & ! [IN]
                          TRACER_MASS(:) ) ! [IN]
    endif

    if ( AD_PREP_sw(I_QTOT) > 0 ) then
!OCL XFILL
       do j  = JSB, JEB
       do i  = ISB, IEB
       do k  = KS, KE
          QTOT(k,i,j) = 1.0_RP - QDRY(k,i,j)
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_QHYD) > 0 ) then
!OCL XFILL
       QHYD(:,:,:) = 0.0_RP
       do iq = QHS, QHE
         QHYD(:,:,:) = QHYD(:,:,:) + QTRC(:,:,:,iq)
       enddo
    endif

    if ( AD_PREP_sw(I_QLIQ) > 0 ) then
!OCL XFILL
       QLIQ(:,:,:) = 0.0_RP
       do iq = QLS, QLE
         QLIQ(:,:,:) = QLIQ(:,:,:) + QTRC(:,:,:,iq)
       enddo
    endif

    if ( AD_PREP_sw(I_QICE) > 0 ) then
!OCL XFILL
       QICE(:,:,:) = 0.0_RP
       do iq = QIS, QIE
         QICE(:,:,:) = QICE(:,:,:) + QTRC(:,:,:,iq)
       enddo
    endif

    if ( AD_PREP_sw(I_LWP) > 0 ) then
       do j  = JSB, JEB
       do i  = ISB, IEB
          LWP(i,j) = 0.0_RP
          do k  = KS, KE
             LWP(i,j) = LWP(i,j) &
                      + QLIQ(k,i,j) * DENS(k,i,j) * ( REAL_FZ(k,i,j)-REAL_FZ(k-1,i,j) ) * 1.E3_RP ! [kg/m2->g/m2]
          enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_IWP) > 0 ) then
       do j  = JSB, JEB
       do i  = ISB, IEB
          IWP(i,j) = 0.0_RP
          do k  = KS, KE
             IWP(i,j) = IWP(i,j) &
                      + QICE(k,i,j) * DENS(k,i,j) * ( REAL_FZ(k,i,j)-REAL_FZ(k-1,i,j) ) * 1.E3_RP ! [kg/m2->g/m2]
          enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_PW) > 0 ) then
       do j  = JSB, JEB
       do i  = ISB, IEB
          PW(i,j) = 0.0_RP
          do k  = KS, KE
             PW(i,j) = PW(i,j) &
                      + QTRC(k,i,j,I_QV) * DENS(k,i,j) * ( REAL_FZ(k,i,j)-REAL_FZ(k-1,i,j) ) * 1.E3_RP ! [kg/m2->g/m2]
          enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_RTOT) > 0 ) then
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
       do j  = JSB, JEB
       do i  = ISB, IEB
       do k  = KS, KE
          CALC_R(RTOT(k,i,j), QDRY(k,i,j), QTRC, k, i, j, iq, Rdry, TRACER_R)
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_CPTOT) > 0 ) then
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
       do j  = JSB, JEB
       do i  = ISB, IEB
       do k  = KS, KE
          CPTOT(k,i,j) = CPdry * QDRY(k,i,j)
          CVTOT(k,i,j) = CVdry * QDRY(k,i,j)
          do iq = 1, QA
             CPTOT(k,i,j) = CPTOT(k,i,j) + QTRC(k,i,j,iq) * TRACER_CP(iq)
             CVTOT(k,i,j) = CVTOT(k,i,j) + QTRC(k,i,j,iq) * TRACER_CV(iq)
          end do
       enddo
       enddo
       enddo

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
       do j  = JSB, JEB
       do i  = ISB, IEB
       do k  = KS, KE
          CPovCV(k,i,j) = CPTOT(k,i,j) / CVTOT(k,i,j)
       enddo
       enddo
       enddo
    endif

!    if ( AD_PREP_sw(I_PRES) > 0 ) then
!       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!       do j = JS, JE
!       do i = IS, IE
!       do k = KS, KE
!          PRES(k,i,j) = P00 * ( RHOT(k,i,j) * RTOT(k,i,j) / P00 )**CPovCV(k,i,j)
!       enddo
!       enddo
!       enddo
!    endif
!
!    if ( AD_PREP_sw(I_TEMP) > 0 ) then
!       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!       do j = JS, JE
!       do i = IS, IE
!       do k = KS, KE
!          TEMP(k,i,j) = PRES(k,i,j) / ( DENS(k,i,j) * RTOT(k,i,j) )
!       enddo
!       enddo
!       enddo
!    endif

    if ( AD_PREP_sw(I_POTL) > 0 ) then
       call HYDROMETEOR_LHV( LHV_local(:,:,:), TEMP(:,:,:) )

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
       do j  = JSB, JEB
       do i  = ISB, IEB
       do k = KS, KE
          POTL(k,i,j) = POTT(k,i,j) &
                      - LHV_local(k,i,j) / CPdry * QLIQ(k,i,j) * POTT(k,i,j) / TEMP(k,i,j)
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_QSAT) > 0 ) then
!       call SATURATION_dens2qsat_all( QSAT(:,:,:), & ! [OUT]
!                                      TEMP(:,:,:), & ! [IN]
!                                      DENS(:,:,:)  ) ! [IN]
    end if

    if ( AD_HIST_id(I_RHA) > 0 ) then
       call SATURATION_psat_all( PSAT(:,:,:), TEMP(:,:,:) )

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
       do j  = JSB, JEB
       do i  = ISB, IEB
       do k  = KS, KE
          RHA(k,i,j) = DENS(k,i,j) * QTRC(k,i,j,I_QV) &
                     / PSAT(k,i,j) * Rvap * TEMP(k,i,j) &
                     * 100.0_RP
       enddo
       enddo
       enddo
    endif

    if ( AD_HIST_id(I_RHL) > 0 ) then
       call SATURATION_psat_liq( PSAT(:,:,:), & ! [OUT]
                                 TEMP(:,:,:)  ) ! [IN]
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
       do j  = JSB, JEB
       do i  = ISB, IEB
       do k  = KS, KE
          RHL(k,i,j) = DENS(k,i,j) * QTRC(k,i,j,I_QV) &
                     / PSAT(k,i,j) * Rvap * TEMP(k,i,j) &
                     * 100.0_RP
       enddo
       enddo
       enddo
    endif

    if ( AD_HIST_id(I_RHI) > 0 ) then
       call SATURATION_psat_ice( PSAT(:,:,:), & ! [OUT]
                                 TEMP(:,:,:)  ) ! [IN]
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
       do j  = JSB, JEB
       do i  = ISB, IEB
       do k  = KS, KE
          RHI(k,i,j) = DENS(k,i,j) * QTRC(k,i,j,I_QV) &
                     / PSAT(k,i,j) * Rvap * TEMP(k,i,j) &
                     * 100.0_RP
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_VOR) > 0 ) then
       ! at x, v, layer
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
       do j = 1, JA-1
       do i = 2, IA
       do k = KS, KE
          UH(k,i,j) = 0.5_RP * ( MOMX(k,i,j)+MOMX(k,i,j+1)+MOMX(k,i-1,j)+MOMX(k,i-1,j+1) ) &
                             / ( DENS(k,i,j)+DENS(k,i,j+1) )
       enddo
       enddo
       enddo

       ! at u, y, layer
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
       do j = 2, JA
       do i = 1, IA-1
       do k = KS, KE
          VH(k,i,j) = 0.5_RP * ( MOMY(k,i,j)+MOMY(k,i+1,j)+MOMY(k,i,j-1)+MOMY(k,i+1,j-1) ) &
                             / ( DENS(k,i,j)+DENS(k,i+1,j) )
       enddo
       enddo
       enddo

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
       do j = 2, JA-1
       do i = 2, IA-1
       do k = KS, KE
          VOR(k,i,j) = ( VH(k,i,j  ) - VH(k,i-1,j  ) ) * RCDX(i) &
                     - ( UH(k,i  ,j) - UH(k,i  ,j-1) ) * RCDY(j)
       enddo
       enddo
       enddo

       !$omp parallel do private(j,k) OMP_SCHEDULE_
       do j = 1, JA
       do k = KS, KE
          VOR(k,1 ,j) = VOR(k,2   ,j)
          VOR(k,IA,j) = VOR(k,IA-1,j)
       enddo
       enddo

       !$omp parallel do private(i,k) OMP_SCHEDULE_
       do i = 1, IA
       do k = KS, KE
          VOR(k,i,1 ) = VOR(k,i,2   )
          VOR(k,i,JA) = VOR(k,i,JA-1)
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_HDIV) > 0 ) then
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
       do j = 2, JA
       do i = 2, IA
       do k = KS, KE
          HDIV(k,i,j) = ( MOMX(k,i,j) - MOMX(k  ,i-1,j  ) ) * RCDX(i) &
                      + ( MOMY(k,i,j) - MOMY(k  ,i  ,j-1) ) * RCDY(j)
       enddo
       enddo
       enddo
       !$omp parallel do private(i,k) OMP_SCHEDULE_
       do i = 1, IA
       do k = KS, KE
          HDIV(k,i,1) = HDIV(k,i,2)
       enddo
       enddo
       !$omp parallel do private(j,k) OMP_SCHEDULE_
       do j = 1, JA
       do k = KS, KE
          HDIV(k,1,j) = HDIV(k,2,j)
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_DIV) > 0 ) then
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          DIV(k,i,j) = ( MOMZ(k,i,j) - MOMZ(k-1,i  ,j  ) ) * ( REAL_FZ(k,i,j)-REAL_FZ(k-1,i,j) ) &
                     + HDIV(k,i,j)
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_Uabs) > 0 ) then
!OCL XFILL
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          Uabs(k,i,j) = sqrt( U(k,i,j) * U(k,i,j) &
                            + V(k,i,j) * V(k,i,j) )
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_DENS_MEAN) > 0 ) then
       call COMM_horizontal_mean( DENS_MEAN(:), DENS(:,:,:) )
    end if

    if ( AD_PREP_sw(I_DENS_PRIM) > 0 ) then
!OCL XFILL
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          DENS_PRIM(k,i,j) = DENS(k,i,j) - DENS_MEAN(k)
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_W_MEAN) > 0 ) then
!OCL XFILL
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          W_PRIM(k,i,j) = W(k,i,j) * DENS(k,i,j)
       enddo
       enddo
       enddo
       call COMM_horizontal_mean( W_MEAN(:), W_PRIM(:,:,:) )
       do k = KS, KE
          W_MEAN(k) = W_MEAN(k) / DENS_MEAN(k)
       enddo
    end if
    if ( AD_PREP_sw(I_W_PRIM) > 0 ) then
!OCL XFILL
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          W_PRIM(k,i,j) = W(k,i,j) - W_MEAN(k)
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_U_MEAN) > 0 ) then
!OCL XFILL
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          U_PRIM(k,i,j) = U(k,i,j) * DENS(k,i,j)
       enddo
       enddo
       enddo
       call COMM_horizontal_mean( U_MEAN(:), U_PRIM(:,:,:) )
       do k = KS, KE
          U_MEAN(k) = U_MEAN(k) / DENS_MEAN(k)
       enddo
    end if
    if ( AD_PREP_sw(I_U_PRIM) > 0 ) then
!OCL XFILL
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          U_PRIM(k,i,j) = U(k,i,j) - U_MEAN(k)
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_V_MEAN) > 0 ) then
!OCL XFILL
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          V_PRIM(k,i,j) = V(k,i,j) * DENS(k,i,j)
       enddo
       enddo
       enddo
       call COMM_horizontal_mean( V_MEAN(:), V_PRIM(:,:,:) )
       do k = KS, KE
          V_MEAN(k) = V_MEAN(k) / DENS_MEAN(k)
       enddo
    end if
    if ( AD_PREP_sw(I_V_PRIM) > 0 ) then
!OCL XFILL
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          V_PRIM(k,i,j) = V(k,i,j) - V_MEAN(k)
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_T_MEAN) > 0 ) then
!OCL XFILL
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          POTT_PRIM(k,i,j) = TEMP(k,i,j) * DENS(k,i,j)
       enddo
       enddo
       enddo
       call COMM_horizontal_mean( T_MEAN(:), POTT_PRIM(:,:,:) )
       do k = KS, KE
          T_MEAN(k) = T_MEAN(k) / DENS_MEAN(k)
       enddo
    end if

    if ( AD_PREP_sw(I_POTT_MEAN) > 0 ) then
       call COMM_horizontal_mean( PT_MEAN(:), RHOT(:,:,:) )
       do k = KS, KE
          PT_MEAN(k) = PT_MEAN(k) / DENS_MEAN(k)
       enddo
    end if
    if ( AD_PREP_sw(I_POTT_PRIM) > 0 ) then
!OCL XFILL
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          POTT_PRIM(k,i,j) = POTT(k,i,j) - PT_MEAN(k)
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_QV_MEAN) > 0 ) then
!OCL XFILL
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          RHOQ(k,i,j) = QTRC(k,i,j,I_QV) * DENS(k,i,j)
       enddo
       enddo
       enddo
       call COMM_horizontal_mean( QV_MEAN(:), RHOQ(:,:,:) )
       do k = KS, KE
          QV_MEAN(k) = QV_MEAN(k) / DENS_MEAN(k)
       enddo
    end if

    if ( AD_PREP_sw(I_QHYD_MEAN) > 0 ) then
!OCL XFILL
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          RHOQ(k,i,j) = QHYD(k,i,j) * DENS(k,i,j)
       enddo
       enddo
       enddo
       call COMM_horizontal_mean( QHYD_MEAN(:), RHOQ(:,:,:) )
       do k = KS, KE
          QHYD_MEAN(k) = QHYD_MEAN(k) / DENS_MEAN(k)
       enddo
    end if

    if ( AD_PREP_sw(I_QLIQ_MEAN) > 0 ) then
!OCL XFILL
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          RHOQ(k,i,j) = QLIQ(k,i,j) * DENS(k,i,j)
       enddo
       enddo
       enddo
       call COMM_horizontal_mean( QLIQ_MEAN(:), RHOQ(:,:,:) )
       do k = KS, KE
          QLIQ_MEAN(k) = QLIQ_MEAN(k) / DENS_MEAN(k)
       enddo
    end if

    if ( AD_PREP_sw(I_QICE_MEAN) > 0 ) then
!OCL XFILL
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          RHOQ(k,i,j) = QICE(k,i,j) * DENS(k,i,j)
       enddo
       enddo
       enddo
       call COMM_horizontal_mean( QICE_MEAN(:), RHOQ(:,:,:) )
       do k = KS, KE
          QICE_MEAN(k) = QICE_MEAN(k) / DENS_MEAN(k)
       enddo
    end if

    if ( AD_PREP_sw(I_W_PRIM2) > 0 ) then
!OCL XFILL
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          W_PRIM2(k,i,j) = W_PRIM(k,i,j) * W_PRIM(k,i,j)
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_PT_W_PRIM) > 0 ) then
!OCL XFILL
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          PT_W_PRIM(k,i,j) = W_PRIM(k,i,j) * POTT_PRIM(k,i,j) * DENS(k,i,j) * CPdry
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_W_PRIM3) > 0 ) then
!OCL XFILL
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          W_PRIM3(k,i,j) = W_PRIM(k,i,j) * W_PRIM(k,i,j) * W_PRIM(k,i,j)
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_TKE_RS) > 0 ) then
!OCL XFILL
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          TKE_RS(k,i,j) = 0.5_RP * ( W_PRIM(k,i,j) * W_PRIM(k,i,j) &
                                   + U_PRIM(k,i,j) * U_PRIM(k,i,j) &
                                   + V_PRIM(k,i,j) * V_PRIM(k,i,j) )
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_ENGP) > 0 ) then
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          ENGP(k,i,j) = DENS(k,i,j) * GRAV * REAL_CZ(k,i,j)
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_ENGK) > 0 ) then
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          ENGK(k,i,j) = 0.5_RP * DENS(k,i,j) * ( W(k,i,j)**2 &
                                               + U(k,i,j)**2 &
                                               + V(k,i,j)**2 )
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_ENGI) > 0 ) then
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
       do j  = JSB, JEB
       do i  = ISB, IEB
       do k = KS, KE
          ENGI(k,i,j) = DENS(k,i,j) * QDRY(k,i,j) * TEMP(k,i,j) * CVdry
       enddo
       enddo
       enddo

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(3)
       do j  = JSB, JEB
       do i  = ISB, IEB
       do k  = KS, KE
       do iq = 1, QA
          ENGI(k,i,j) = ENGI(k,i,j) &
                      + DENS(k,i,j) * QTRC(k,i,j,iq) * TEMP(k,i,j) * TRACER_CV(iq)
       enddo
       enddo
       enddo
       enddo

       if ( I_QV > 0 ) then
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(3)
       do j  = JSB, JEB
       do i  = ISB, IEB
       do k  = KS, KE
          ENGI(k,i,j) = ENGI(k,i,j) &
                      + DENS(k,i,j) * QTRC(k,i,j,I_QV) * LHV ! Latent Heat [vapor->liquid]
       enddo
       enddo
       enddo
       end if

       do iq = QIS, QIE
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(3)
       do j  = JSB, JEB
       do i  = ISB, IEB
       do k  = KS, KE
          ENGI(k,i,j) = ENGI(k,i,j) &
                      - DENS(k,i,j) * QTRC(k,i,j,iq) * LHF ! Latent Heat [ice->liquid]
       enddo
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_ENGT) > 0 ) then
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
       do j  = JSB, JEB
       do i  = ISB, IEB
       do k  = KS, KE
          ENGT(k,i,j) = ENGP(k,i,j) + ENGK(k,i,j) + ENGI(k,i,j)
       enddo
       enddo
       enddo
    endif

    call HIST_in( W    (:,:,:), 'W',     'velocity w',             'm/s'    )
    call HIST_in( U    (:,:,:), 'U',     'velocity u',             'm/s'    )
    call HIST_in( V    (:,:,:), 'V',     'velocity v',             'm/s'    )
    call HIST_in( POTT (:,:,:), 'PT',    'potential temp.',        'K'      )

    call HIST_in( QDRY (:,:,:), 'QDRY',  'dry air',                'kg/kg'  )
    call HIST_in( QTOT (:,:,:), 'QTOT',  'total water',            'kg/kg'  )
    call HIST_in( QHYD (:,:,:), 'QHYD',  'total hydrometeors',     'kg/kg'  )
    call HIST_in( QLIQ (:,:,:), 'QLIQ',  'total liquid water',     'kg/kg'  )
    call HIST_in( QICE (:,:,:), 'QICE',  'total ice water',        'kg/kg'  )

    call HIST_in( LWP  (:,:),   'LWP',   'liquid water path',      'g/m2'   )
    call HIST_in( IWP  (:,:),   'IWP',   'ice    water path',      'g/m2'   )
    call HIST_in( PW   (:,:),   'PW',    'precipitable water',     'g/m2'   )

    call HIST_in( RTOT (:,:,:), 'RTOT',  'Total gas constant',     'J/kg/K' )
    call HIST_in( CPTOT(:,:,:), 'CPTOT', 'Total heat capacity',    'J/kg/K' )
    call HIST_in( PRES (:,:,:), 'PRES',  'pressure',               'Pa'     )
    call HIST_in( TEMP (:,:,:), 'T',     'temperature',            'K'      )

    call HIST_in( POTL (:,:,:), 'LWPT',  'liq. potential temp.',   'K'      )
    call HIST_in( RHA  (:,:,:), 'RHA',   'relative humidity(liq+ice)','%'   )
    call HIST_in( RHL  (:,:,:), 'RH' ,   'relative humidity(liq)', '%'      )
    call HIST_in( RHI  (:,:,:), 'RHI',   'relative humidity(ice)', '%'      )

    call HIST_in( VOR  (:,:,:), 'VOR',   'vertical vorticity',     '1/s'    )
    call HIST_in( DIV  (:,:,:), 'DIV',   'divergence',             '1/s'    )
    call HIST_in( HDIV (:,:,:), 'HDIV',  'horizontal divergence',  '1/s'    )
    call HIST_in( Uabs (:,:,:), 'Uabs',  'absolute velocity',      'm/s'    )

    call HIST_in( DENS_MEAN(:), 'DENS_MEAN', 'horiz. mean of density',    'kg/m3' )
    call HIST_in( W_MEAN   (:), 'W_MEAN',    'horiz. mean of w',          'm/s' )
    call HIST_in( U_MEAN   (:), 'U_MEAN',    'horiz. mean of u',          'm/s' )
    call HIST_in( V_MEAN   (:), 'V_MEAN',    'horiz. mean of v',          'm/s' )
    call HIST_in( PT_MEAN  (:), 'PT_MEAN',   'horiz. mean of pot.',       'K' )
    call HIST_in( T_MEAN   (:), 'T_MEAN',    'horiz. mean of t',          'K' )
    call HIST_in( QV_MEAN  (:), 'QV_MEAN',   'horiz. mean of QV',         '1' )
    call HIST_in( QHYD_MEAN(:), 'QHYD_MEAN', 'horiz. mean of QHYD',       '1' )
    call HIST_in( QLIQ_MEAN(:), 'QLIQ_MEAN', 'horiz. mean of QLIQ',       '1' )
    call HIST_in( QICE_MEAN(:), 'QICE_MEAN', 'horiz. mean of QICE',       '1' )

    call HIST_in( DENS_PRIM(:,:,:), 'DENS_PRIM', 'horiz. deviation of density',    'kg/m3' )
    call HIST_in( W_PRIM   (:,:,:), 'W_PRIM',    'horiz. deviation of w',          'm/s'   )
    call HIST_in( U_PRIM   (:,:,:), 'U_PRIM',    'horiz. deviation of u',          'm/s'   )
    call HIST_in( V_PRIM   (:,:,:), 'V_PRIM',    'horiz. deviation of v',          'm/s'   )
    call HIST_in( POTT_PRIM(:,:,:), 'PT_PRIM',   'horiz. deviation of pot. temp.', 'K'     )
    call HIST_in( W_PRIM2  (:,:,:), 'W_PRIM2',   'variance of w',                  'm2/s2' )
    call HIST_in( PT_W_PRIM(:,:,:), 'PT_W_PRIM', 'resolved scale heat flux',       'W/s'   )
    call HIST_in( W_PRIM3  (:,:,:), 'W_PRIM3',   'skewness of w',                  'm3/s3' )
    call HIST_in( TKE_RS   (:,:,:), 'TKE_RS',    'resolved scale TKE',             'm2/s2' )

    call HIST_in( ENGT (:,:,:), 'ENGT',  'total energy',           'J/m3'   )
    call HIST_in( ENGP (:,:,:), 'ENGP',  'potential energy',       'J/m3'   )
    call HIST_in( ENGK (:,:,:), 'ENGK',  'kinetic energy',         'J/m3'   )
    call HIST_in( ENGI (:,:,:), 'ENGI',  'internal energy',        'J/m3'   )

    if (      AD_PREP_sw(I_CAPE) > 0 &
         .OR. AD_PREP_sw(I_CIN)  > 0 &
         .OR. AD_PREP_sw(I_LCL)  > 0 &
         .OR. AD_PREP_sw(I_LFC)  > 0 &
         .OR. AD_PREP_sw(I_LNB)  > 0 ) then

       call ADIABAT_cape( KS,               & ! [IN]
                          DENS   (:,:,:),   & ! [IN]
                          TEMP   (:,:,:),   & ! [IN]
                          PRES   (:,:,:),   & ! [IN]
                          QTRC   (:,:,:,:), & ! [IN]
                          REAL_CZ(:,:,:),   & ! [IN]
                          REAL_FZ(:,:,:),   & ! [IN]
                          CAPE   (:,:),     & ! [OUT]
                          CIN    (:,:),     & ! [OUT]
                          LCL    (:,:),     & ! [OUT]
                          LFC    (:,:),     & ! [OUT]
                          LNB    (:,:)      ) ! [OUT]

    endif

    call HIST_in( CAPE(:,:), 'CAPE', 'convection avail. pot. energy', 'm2/s2' )
    call HIST_in( CIN (:,:), 'CIN',  'convection inhibition',         'm2/s2' )
    call HIST_in( LCL (:,:), 'LCL',  'lifted condensation level',     'm'     )
    call HIST_in( LFC (:,:), 'LFC',  'level of free convection',      'm'     )
    call HIST_in( LNB (:,:), 'LNB',  'level of neutral buoyancy',     'm'     )

    if ( AD_PREP_sw(I_PBLH) > 0 ) then
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          fact = 1.0_RP
          if ( I_QV > 0 ) fact = fact + 0.61_RP * QTRC(k,i,j,I_QV)
          if ( I_QC > 0 ) fact = fact - 1.61 * QTRC(k,i,j,I_QC)
          POTTv(k,i,j) = POTT(k,i,j) * fact
       enddo
       enddo
       enddo

       do j = JS, JE
       do i = IS, IE
          PBLH(i,j) = REAL_CZ(KS,i,j) - REAL_FZ(KS-1,i,j)

          do k = KS+1, KE
             if ( POTTv(k,i,j) > POTTv(KS,i,j) ) then
                fact = ( POTTv(KS,i,j) - POTTv(k-1,i,j) ) &
                     / ( POTTv(k,i,j)  - POTTv(k-1,i,j) )

                PBLH(i,j) = REAL_CZ(k-1,i,j) - REAL_FZ(KS-1,i,j) &
                          + fact * ( REAL_CZ(k,i,j) - REAL_CZ(k-1,i,j) )

                exit
             endif
          enddo
       enddo
       enddo
    endif
    call HIST_in( PBLH(:,:), 'PBLH', 'PBL height', 'm' )

    if ( AD_PREP_sw(I_MSE) > 0 ) then
       call HYDROMETEOR_LHV( LHV_local(:,:,:), TEMP(:,:,:) )

       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          MSE(k,i,j) = CPTOT(k,i,j) * TEMP(k,i,j)                    &
                     + GRAV * ( REAL_CZ(k,i,j) - REAL_FZ(KS-1,i,j) ) &
                     + LHV_local(k,i,j) * QTRC(k,i,j,I_QV)
       enddo
       enddo
       enddo
    endif
    call HIST_in( MSE(:,:,:), 'MSE', 'moist static energy', 'm2/s2' )

    do j = JS, JE
    do i = IS, IE
       PREC(i,j) = SFLX_rain_CP(i,j)                     &
                 + SFLX_rain_MP(i,j) + SFLX_snow_MP(i,j)
    enddo
    enddo
    call HIST_in( PREC(:,:), 'PREC', 'surface precipitation rate (total)', 'kg/m2/s' )

    return
  end subroutine ATMOS_vars_history

  !-----------------------------------------------------------------------------
  !> Budget monitor for atmosphere
  subroutine ATMOS_vars_total
    use scale_const, only: &
       GRAV  => CONST_GRAV,  &
       CVdry => CONST_CVdry
    use scale_grid_real, only: &
       REAL_CZ
    use scale_rm_statistics, only: &
       STATISTICS_checktotal, &
       STAT_total
    use scale_atmos_thermodyn, only: &
       THERMODYN_qd        => ATMOS_THERMODYN_qd,        &
       THERMODYN_temp_pres => ATMOS_THERMODYN_temp_pres
    use scale_atmos_hydrometeor, only: &
       I_QV, &
       QIS,  &
       QIE,  &
       LHV,  &
       LHF
    implicit none

    real(RP) :: W   (KA,IA,JA) ! velocity w at cell center [m/s]
    real(RP) :: U   (KA,IA,JA) ! velocity u at cell center [m/s]
    real(RP) :: V   (KA,IA,JA) ! velocity v at cell center [m/s]

    real(RP) :: QDRY(KA,IA,JA) ! dry air     [kg/kg]
    real(RP) :: PRES(KA,IA,JA) ! pressure    [Pa]
    real(RP) :: TEMP(KA,IA,JA) ! temperature [K]

    real(RP) :: ENGT(KA,IA,JA) ! total     energy [J/m3]
    real(RP) :: ENGP(KA,IA,JA) ! potential energy [J/m3]
    real(RP) :: ENGK(KA,IA,JA) ! kinetic   energy [J/m3]
    real(RP) :: ENGI(KA,IA,JA) ! internal  energy [J/m3]

    real(RP) :: RHOQ(KA,IA,JA)

    real(RP) :: total ! dummy
    integer :: i, j, k, iq
    !---------------------------------------------------------------------------

    if ( STATISTICS_checktotal ) then

       call STAT_total( total, DENS(:,:,:), VAR_NAME(I_DENS) )
       call STAT_total( total, MOMZ(:,:,:), VAR_NAME(I_MOMZ) )
       call STAT_total( total, MOMX(:,:,:), VAR_NAME(I_MOMX) )
       call STAT_total( total, MOMY(:,:,:), VAR_NAME(I_MOMY) )
       call STAT_total( total, RHOT(:,:,:), VAR_NAME(I_RHOT) )

       do iq = 1, QA
          RHOQ(:,:,:) = DENS(:,:,:) * QTRC(:,:,:,iq)

          call STAT_total( total, RHOQ(:,:,:), TRACER_NAME(iq) )
       enddo

       call THERMODYN_qd( QDRY(:,:,:),   & ! [OUT]
                          QTRC(:,:,:,:), & ! [IN]
                          TRACER_MASS(:) ) ! [IN]

       call THERMODYN_temp_pres( TEMP(:,:,:),   & ! [OUT]
                                 PRES(:,:,:),   & ! [OUT]
                                 DENS(:,:,:),   & ! [IN]
                                 RHOT(:,:,:),   & ! [IN]
                                 QTRC(:,:,:,:), & ! [IN]
                                 TRACER_CV(:),  & ! [IN]
                                 TRACER_R(:),   & ! [IN]
                                 TRACER_MASS(:) ) ! [IN]

       RHOQ(KS:KE,IS:IE,JS:JE) = DENS(KS:KE,IS:IE,JS:JE) * QDRY (KS:KE,IS:IE,JS:JE)

       call STAT_total( total, RHOQ(:,:,:), 'QDRY' )

       RHOQ(KS:KE,IS:IE,JS:JE) = DENS(KS:KE,IS:IE,JS:JE) * ( 1.0_RP - QDRY (KS:KE,IS:IE,JS:JE) ) ! Qtotal

       call STAT_total( total, RHOQ(:,:,:), 'QTOT' )

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          W(k,i,j) = 0.5_RP * ( MOMZ(k-1,i,j)+MOMZ(k,i,j) ) / DENS(k,i,j)
          U(k,i,j) = 0.5_RP * ( MOMX(k,i-1,j)+MOMX(k,i,j) ) / DENS(k,i,j)
          V(k,i,j) = 0.5_RP * ( MOMY(k,i,j-1)+MOMY(k,i,j) ) / DENS(k,i,j)

          ENGP(k,i,j) = DENS(k,i,j) * GRAV * REAL_CZ(k,i,j)

          ENGK(k,i,j) = 0.5_RP * DENS(k,i,j) * ( W(k,i,j)**2 &
                                               + U(k,i,j)**2 &
                                               + V(k,i,j)**2 )

          ENGI(k,i,j) = DENS(k,i,j) * QDRY(k,i,j) * TEMP(k,i,j) * CVdry
          do iq = 1, QA
             ENGI(k,i,j) = ENGI(k,i,j) &
                         + DENS(k,i,j) * QTRC(k,i,j,iq) * TEMP(k,i,j) * TRACER_CV(iq)
          enddo

          if ( I_QV > 0 ) then
             ENGI(k,i,j) = ENGI(k,i,j) + DENS(k,i,j) * QTRC(k,i,j,I_QV) * LHV ! Latent Heat [vapor->liquid]
          end if

          do iq = QIS, QIE
             ENGI(k,i,j) = ENGI(k,i,j) - DENS(k,i,j) * QTRC(k,i,j,iq) * LHF ! Latent Heat [ice->liquid]
          enddo

          ENGT(k,i,j) = ENGP(k,i,j) + ENGK(k,i,j) + ENGI(k,i,j)
       enddo
       enddo
       enddo

       call STAT_total( total, ENGP(:,:,:), 'ENGP' )
       call STAT_total( total, ENGK(:,:,:), 'ENGK' )
       call STAT_total( total, ENGI(:,:,:), 'ENGI' )
       call STAT_total( total, ENGT(:,:,:), 'ENGT' )

    endif

    return
  end subroutine ATMOS_vars_total

  !-----------------------------------------------------------------------------
  !> Calc diagnostic variables
  subroutine ATMOS_vars_diagnostics
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    use scale_atmos_thermodyn, only: &
       THERMODYN_temp_pres => ATMOS_THERMODYN_temp_pres
    implicit none

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Calc diagnostics'

    call THERMODYN_temp_pres( TEMP(:,:,:),   & ! [OUT]
                              PRES(:,:,:),   & ! [OUT]
                              DENS(:,:,:),   & ! [IN]
                              RHOT(:,:,:),   & ! [IN]
                              QTRC(:,:,:,:), & ! [IN]
                              TRACER_CV(:),  & ! [IN]
                              TRACER_R(:),   & ! [IN]
                              TRACER_MASS(:) ) ! [IN]

!OCL XFILL
    do j = 1, JA
    do i = 1, IA
    do k = KS+1, KE-1
       W(k,i,j) = 0.5_RP * ( MOMZ(k-1,i,j)+MOMZ(k,i,j) ) / DENS(k,i,j)
    enddo
    enddo
    enddo
!OCL XFILL
    do j = 1, JA
    do i = 1, IA
       W(KS,i,j) = 0.5_RP * (               MOMZ(KS,i,j) ) / DENS(KS,i,j)
    enddo
    enddo
!OCL XFILL
    do j = 1, JA
    do i = 1, IA
       W(KE,i,j) = 0.5_RP * ( MOMZ(KE-1,i,j)             ) / DENS(KE,i,j)
    enddo
    enddo

!OCL XFILL
    do j = 1, JA
    do i = 2, IA
    do k = KS, KE
       U(k,i,j) = 0.5_RP * ( MOMX(k,i-1,j)+MOMX(k,i,j) ) / DENS(k,i,j)
    enddo
    enddo
    enddo
!OCL XFILL
    do j = 1, JA
    do k = KS, KE
       U(k,1,j) = MOMX(k,1,j) / DENS(k,1,j)
    enddo
    enddo

!OCL XFILL
    do j = 2, JA
    do i = 1, IA
    do k = KS, KE
       V(k,i,j) = 0.5_RP * ( MOMY(k,i,j-1)+MOMY(k,i,j) ) / DENS(k,i,j)
    enddo
    enddo
    enddo
!OCL XFILL
    do i = 1, IA
    do k = KS, KE
       V(k,i,1) = MOMY(k,i,1) / DENS(k,i,1)
    enddo
    enddo

    !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
    do j  = 1, JA
    do i  = 1, IA
       W(   1:KS-1,i,j) = W(KS,i,j)
       U(   1:KS-1,i,j) = U(KS,i,j)
       V(   1:KS-1,i,j) = V(KS,i,j)
       W(KE+1:KA,  i,j) = W(KE,i,j)
       U(KE+1:KA,  i,j) = U(KE,i,j)
       V(KE+1:KA,  i,j) = V(KE,i,j)
    enddo
    enddo

    call COMM_vars8( U(:,:,:), 1 )
    call COMM_vars8( V(:,:,:), 2 )
    call COMM_wait ( U(:,:,:), 1, .false. )
    call COMM_wait ( V(:,:,:), 2, .false. )

!OCL XFILL
    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       POTT(k,i,j) = RHOT(k,i,j) / DENS(k,i,j)
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_vars_diagnostics

  !-----------------------------------------------------------------------------
  !> monitor output
  subroutine ATMOS_vars_monitor
    use scale_process, only: &
       PRC_myrank
    use scale_const, only: &
       GRAV  => CONST_GRAV,  &
       CVdry => CONST_CVdry
    use scale_grid, only: &
       RFDX => GRID_RFDX, &
       RFDY => GRID_RFDY
    use scale_grid_real, only: &
       REAL_CZ, &
       REAL_FZ
    use scale_gridtrans, only: &
       MAPF => GTRANS_MAPF, &
       I_UY, &
       I_XV
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    use scale_rm_statistics, only: &
       STATISTICS_checktotal, &
       STAT_total,            &
       STAT_detail
    use scale_monitor, only: &
       MONIT_put, &
       MONIT_in
    use scale_time, only: &
       TIME_DTSEC_ATMOS_DYN
    use scale_atmos_thermodyn, only: &
       THERMODYN_qd        => ATMOS_THERMODYN_qd,        &
       THERMODYN_temp_pres => ATMOS_THERMODYN_temp_pres
    use scale_atmos_hydrometeor, only: &
       I_QV, &
       QIS,  &
       QIE,  &
       LHV,  &
       LHF
    use mod_atmos_phy_mp_vars, only: &
       SFLX_rain => ATMOS_PHY_MP_SFLX_rain, &
       SFLX_snow => ATMOS_PHY_MP_SFLX_snow
    use mod_atmos_phy_rd_vars, only: &
       SFLX_LW_up   => ATMOS_PHY_RD_SFLX_LW_up,   &
       SFLX_LW_dn   => ATMOS_PHY_RD_SFLX_LW_dn,   &
       SFLX_SW_up   => ATMOS_PHY_RD_SFLX_SW_up,   &
       SFLX_SW_dn   => ATMOS_PHY_RD_SFLX_SW_dn,   &
       TOAFLX_LW_up => ATMOS_PHY_RD_TOAFLX_LW_up, &
       TOAFLX_LW_dn => ATMOS_PHY_RD_TOAFLX_LW_dn, &
       TOAFLX_SW_up => ATMOS_PHY_RD_TOAFLX_SW_up, &
       TOAFLX_SW_dn => ATMOS_PHY_RD_TOAFLX_SW_dn
    use mod_atmos_phy_sf_vars, only: &
       SFLX_SH   => ATMOS_PHY_SF_SFLX_SH, &
       SFLX_LH   => ATMOS_PHY_SF_SFLX_LH, &
       SFLX_QTRC => ATMOS_PHY_SF_SFLX_QTRC
    implicit none

    real(RP) :: QDRY(KA,IA,JA) ! dry air         [kg/kg]
    real(RP) :: RHOQ(KA,IA,JA) ! DENS * tracer   [kg/m3]
    real(RP) :: PRCP(IA,JA)    ! rain + snow     [kg/m2/s]

    real(RP) :: ENGT(KA,IA,JA) ! total     energy [J/m3]
    real(RP) :: ENGP(KA,IA,JA) ! potential energy [J/m3]
    real(RP) :: ENGK(KA,IA,JA) ! kinetic   energy [J/m3]
    real(RP) :: ENGI(KA,IA,JA) ! internal  energy [J/m3]

    real(RP) :: ENGFLXT      (IA,JA) ! total flux             [J/m2/s]
    real(RP) :: SFLX_RD_net  (IA,JA) ! net SFC radiation flux [J/m2/s]
    real(RP) :: TOAFLX_RD_net(IA,JA) ! net TOA radiation flux [J/m2/s]

    real(RP)               :: WORK (KA,IA,JA,3)
    character(len=H_SHORT) :: WNAME(3)
    real(RP)               :: CFLMAX

    integer  :: k, i, j, iq
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Monitor'

    call MONIT_in( DENS(:,:,:), VAR_NAME(I_DENS), VAR_DESC(I_DENS), VAR_UNIT(I_DENS), ndim=3, isflux=.false. )
    call MONIT_in( MOMZ(:,:,:), VAR_NAME(I_MOMZ), VAR_DESC(I_MOMZ), VAR_UNIT(I_MOMZ), ndim=3, isflux=.false. )
    call MONIT_in( MOMX(:,:,:), VAR_NAME(I_MOMX), VAR_DESC(I_MOMX), VAR_UNIT(I_MOMX), ndim=3, isflux=.false. )
    call MONIT_in( MOMY(:,:,:), VAR_NAME(I_MOMY), VAR_DESC(I_MOMY), VAR_UNIT(I_MOMY), ndim=3, isflux=.false. )
    call MONIT_in( RHOT(:,:,:), VAR_NAME(I_RHOT), VAR_DESC(I_RHOT), VAR_UNIT(I_RHOT), ndim=3, isflux=.false. )

    !##### Mass Budget #####

    do iq = 1, QA
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          RHOQ(k,i,j) = DENS(k,i,j) * QTRC(k,i,j,iq)
       enddo
       enddo
       enddo

       call MONIT_in( RHOQ(:,:,:), TRACER_NAME(iq), TRACER_DESC(iq), TRACER_UNIT(iq), ndim=3, isflux=.false. )
    enddo

    ! total dry airmass

    call THERMODYN_qd( QDRY(:,:,:),   & ! [OUT]
                       QTRC(:,:,:,:), & ! [IN]
                       TRACER_MASS(:) ) ! [IN]

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       RHOQ(k,i,j) = DENS(k,i,j) * QDRY (k,i,j)
    enddo
    enddo
    enddo
    call MONIT_put( AD_MONIT_id(I_QDRY), RHOQ(:,:,:) )

    ! total vapor,liquid,solid tracers
    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       RHOQ(k,i,j) = DENS(k,i,j) * ( 1.0_RP - QDRY (k,i,j) )
    enddo
    enddo
    enddo
    call MONIT_put( AD_MONIT_id(I_QTOT), RHOQ(:,:,:) )

    ! total evapolation
    call MONIT_put( AD_MONIT_id(I_EVAP), SFLX_QTRC(:,:,I_QV) )

    ! total precipitation
!OCL XFILL
    do j = JS, JE
    do i = IS, IE
       PRCP(i,j) = SFLX_rain(i,j) + SFLX_snow(i,j)
    enddo
    enddo
    call MONIT_put( AD_MONIT_id(I_PRCP), PRCP(:,:) )

    !##### Energy Budget #####

    call THERMODYN_temp_pres( TEMP(:,:,:),   & ! [OUT]
                              PRES(:,:,:),   & ! [OUT]
                              DENS(:,:,:),   & ! [IN]
                              RHOT(:,:,:),   & ! [IN]
                              QTRC(:,:,:,:), & ! [IN]
                              TRACER_CV(:),  & ! [IN]
                              TRACER_R(:),   & ! [IN]
                              TRACER_MASS(:) ) ! [IN]

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       ENGP(k,i,j) = DENS(k,i,j) * GRAV * REAL_CZ(k,i,j)

       ENGK(k,i,j) = 0.5_RP * DENS(k,i,j) * ( W(k,i,j)**2 &
                                            + U(k,i,j)**2 &
                                            + V(k,i,j)**2 )

       ENGI(k,i,j) = DENS(k,i,j) * QDRY(k,i,j) * TEMP(k,i,j) * CVdry
       do iq = 1, QA
          ENGI(k,i,j) = ENGI(k,i,j) &
                      + DENS(k,i,j) * QTRC(k,i,j,iq) * TEMP(k,i,j) * TRACER_CV(iq)
       enddo

       if ( I_QV > 0 ) then
          ENGI(k,i,j) = ENGI(k,i,j) + DENS(k,i,j) * QTRC(k,i,j,I_QV) * LHV ! Latent Heat [vapor->liquid]
       end if

       do iq = QIS, QIE
          ENGI(k,i,j) = ENGI(k,i,j) - DENS(k,i,j) * QTRC(k,i,j,iq) * LHF ! Latent Heat [ice->liquid]
       enddo
    enddo
    enddo
    enddo

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       ENGT(k,i,j) = ENGP(k,i,j) + ENGK(k,i,j) + ENGI(k,i,j)
    enddo
    enddo
    enddo

!OCL XFILL
    do j = JS, JE
    do i = IS, IE
       SFLX_RD_net(i,j) = ( SFLX_LW_up(i,j) - SFLX_LW_dn(i,j) ) &
                        + ( SFLX_SW_up(i,j) - SFLX_SW_dn(i,j) )

       TOAFLX_RD_net(i,j) = ( TOAFLX_LW_up(i,j) - TOAFLX_LW_dn(i,j) ) &
                          + ( TOAFLX_SW_up(i,j) - TOAFLX_SW_dn(i,j) )
    enddo
    enddo

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
    do j = JS, JE
    do i = IS, IE
       ENGFLXT(i,j) = SFLX_SH(i,j) + SFLX_LH(i,j) &
                    + SFLX_RD_net(i,j) - TOAFLX_RD_net(i,j)
    enddo
    enddo

    call MONIT_put( AD_MONIT_id(I_ENGP), ENGP(:,:,:) )
    call MONIT_put( AD_MONIT_id(I_ENGK), ENGK(:,:,:) )
    call MONIT_put( AD_MONIT_id(I_ENGI), ENGI(:,:,:) )
    call MONIT_put( AD_MONIT_id(I_ENGT), ENGT(:,:,:) )

    call MONIT_put( AD_MONIT_id(I_ENGFLXT), ENGFLXT(:,:) )


    call MONIT_put( AD_MONIT_id(I_ENGSFC_SH), SFLX_SH      (:,:) )
    call MONIT_put( AD_MONIT_id(I_ENGSFC_LH), SFLX_LH      (:,:) )
    call MONIT_put( AD_MONIT_id(I_ENGSFC_RD), SFLX_RD_net  (:,:) )
    call MONIT_put( AD_MONIT_id(I_ENGTOA_RD), TOAFLX_RD_net(:,:) )

    call MONIT_put( AD_MONIT_id(I_ENGSFC_LW_up), SFLX_LW_up(:,:) )
    call MONIT_put( AD_MONIT_id(I_ENGSFC_LW_dn), SFLX_LW_dn(:,:) )
    call MONIT_put( AD_MONIT_id(I_ENGSFC_SW_up), SFLX_SW_up(:,:) )
    call MONIT_put( AD_MONIT_id(I_ENGSFC_SW_dn), SFLX_SW_dn(:,:) )

    call MONIT_put( AD_MONIT_id(I_ENGTOA_LW_up), TOAFLX_LW_up(:,:) )
    call MONIT_put( AD_MONIT_id(I_ENGTOA_LW_dn), TOAFLX_LW_dn(:,:) )
    call MONIT_put( AD_MONIT_id(I_ENGTOA_SW_up), TOAFLX_SW_up(:,:) )
    call MONIT_put( AD_MONIT_id(I_ENGTOA_SW_dn), TOAFLX_SW_dn(:,:) )

    if ( ATMOS_VARS_CHECKRANGE ) then
!OCL XFILL
       WORK(:,:,:,1) = W(:,:,:)
!OCL XFILL
       WORK(:,:,:,2) = U(:,:,:)
!OCL XFILL
       WORK(:,:,:,3) = V(:,:,:)

       WNAME(1) = "W"
       WNAME(2) = "U"
       WNAME(3) = "V"

       call STAT_detail( WORK(:,:,:,:), WNAME(:) )
    endif

    if ( ATMOS_VARS_CHECKCFL > 0.0_RP ) then
!OCL XFILL
       WORK(:,:,:,:) = 0.0_RP

       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          WORK(k,i,j,1) = 0.5_RP * abs(MOMZ(k,i,j)) / ( DENS(k+1,i,j) + DENS(k,i,j) ) &
                        * TIME_DTSEC_ATMOS_DYN / ( REAL_CZ(k+1,i,j) - REAL_CZ(k,i,j) )
          WORK(k,i,j,2) = 0.5_RP * abs(MOMX(k,i,j)) / ( DENS(k,i+1,j) + DENS(k,i,j) ) &
                        * TIME_DTSEC_ATMOS_DYN * RFDX(i) * MAPF(i,j,1,I_UY)
          WORK(k,i,j,3) = 0.5_RP * abs(MOMY(k,i,j)) / ( DENS(k,i,j+1) + DENS(k,i,j) ) &
                        * TIME_DTSEC_ATMOS_DYN * RFDY(j) * MAPF(i,j,2,I_XV)
       enddo
       enddo
       enddo

       CFLMAX = maxval( WORK(:,:,:,:) )
       if ( CFLMAX > ATMOS_VARS_CHECKCFL ) then
          if( IO_L ) write(IO_FID_LOG,*) "*** [ATMOS_vars_monitor] Courant number exceeded the upper limit. : ", CFLMAX
                     write(*,*)          "*** [ATMOS_vars_monitor] Courant number exceeded the upper limit. : ", CFLMAX, &
                                         ", rank = ", PRC_myrank

          WNAME(1) = "Courant num. Z"
          WNAME(2) = "Courant num. X"
          WNAME(3) = "Courant num. Y"

          call STAT_detail( WORK(:,:,:,:), WNAME(:), supress_globalcomm=.true. )
       endif
    endif

    return
  end subroutine ATMOS_vars_monitor

  !-----------------------------------------------------------------------------
  !> Create atmospheric restart file
  subroutine ATMOS_vars_restart_create
    use scale_time, only: &
       TIME_gettimelabel
    use scale_fileio, only: &
       FILEIO_create
    use mod_atmos_admin, only: &
       ATMOS_sw_dyn,      &
       ATMOS_sw_phy_mp,   &
       ATMOS_sw_phy_ae,   &
       ATMOS_sw_phy_ch,   &
       ATMOS_sw_phy_rd,   &
       ATMOS_sw_phy_sf,   &
       ATMOS_sw_phy_tb,   &
       ATMOS_sw_phy_cp
    use mod_atmos_dyn_vars, only: &
       ATMOS_DYN_vars_restart_create
    use mod_atmos_phy_mp_vars, only: &
       ATMOS_PHY_MP_vars_restart_create
    use mod_atmos_phy_ae_vars, only: &
       ATMOS_PHY_AE_vars_restart_create
    use mod_atmos_phy_ch_vars, only: &
       ATMOS_PHY_CH_vars_restart_create
    use mod_atmos_phy_rd_vars, only: &
       ATMOS_PHY_RD_vars_restart_create
    use mod_atmos_phy_sf_vars, only: &
       ATMOS_PHY_SF_vars_restart_create
    use mod_atmos_phy_tb_vars, only: &
       ATMOS_PHY_TB_vars_restart_create
    use mod_atmos_phy_cp_vars, only: &
       ATMOS_PHY_CP_vars_restart_create
#ifdef _SDM
    use scale_atmos_phy_mp_sdm, only: &
       sd_rest_flg_out, &
       ATMOS_PHY_MP_sdm_restart_create
    use scale_time, only: &
       NOWSEC => TIME_NOWSEC
#endif
    implicit none

    character(len=20)     :: timelabel
    character(len=H_LONG) :: basename

    integer :: iq
    !---------------------------------------------------------------------------

#ifdef _SDM
    if( sd_rest_flg_out ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Output random number for SDM ***'
       call ATMOS_PHY_MP_sdm_restart_create(NOWSEC)
    endif
#endif

    if ( ATMOS_RESTART_OUT_BASENAME /= '' ) then

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Output restart file (ATMOS) ***'

       if ( ATMOS_RESTART_OUT_POSTFIX_TIMELABEL ) then
          call TIME_gettimelabel( timelabel )
          basename = trim(ATMOS_RESTART_OUT_BASENAME)//'_'//trim(timelabel)
       else
          basename = trim(ATMOS_RESTART_OUT_BASENAME)
       endif

       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(basename)

       call FILEIO_create( restart_fid,                                               & ! [OUT]
                           basename, ATMOS_RESTART_OUT_TITLE, ATMOS_RESTART_OUT_DTYPE ) ! [IN]

       allocate( VAR_ID(VMAX+QA) )
    endif

    if( ATMOS_sw_dyn )    call ATMOS_DYN_vars_restart_create
    if( ATMOS_sw_phy_mp ) call ATMOS_PHY_MP_vars_restart_create
    if( ATMOS_sw_phy_ae ) call ATMOS_PHY_AE_vars_restart_create
    if( ATMOS_sw_phy_ch ) call ATMOS_PHY_CH_vars_restart_create
    if( ATMOS_sw_phy_rd ) call ATMOS_PHY_RD_vars_restart_create
    if( ATMOS_sw_phy_sf ) call ATMOS_PHY_SF_vars_restart_create
    if( ATMOS_sw_phy_tb ) call ATMOS_PHY_TB_vars_restart_create
    if( ATMOS_sw_phy_cp ) call ATMOS_PHY_CP_vars_restart_create

    return
  end subroutine ATMOS_vars_restart_create

  !-----------------------------------------------------------------------------
  !> Exit netCDF define mode
  subroutine ATMOS_vars_restart_enddef
    use scale_fileio, only: &
       FILEIO_enddef
    use mod_atmos_admin, only: &
       ATMOS_sw_dyn,      &
       ATMOS_sw_phy_mp,   &
       ATMOS_sw_phy_ae,   &
       ATMOS_sw_phy_ch,   &
       ATMOS_sw_phy_rd,   &
       ATMOS_sw_phy_sf,   &
       ATMOS_sw_phy_tb,   &
       ATMOS_sw_phy_cp
    use mod_atmos_dyn_vars, only: &
       ATMOS_DYN_vars_restart_enddef
    use mod_atmos_phy_mp_vars, only: &
       ATMOS_PHY_MP_vars_restart_enddef
    use mod_atmos_phy_ae_vars, only: &
       ATMOS_PHY_AE_vars_restart_enddef
    use mod_atmos_phy_ch_vars, only: &
       ATMOS_PHY_CH_vars_restart_enddef
    use mod_atmos_phy_rd_vars, only: &
       ATMOS_PHY_RD_vars_restart_enddef
    use mod_atmos_phy_sf_vars, only: &
       ATMOS_PHY_SF_vars_restart_enddef
    use mod_atmos_phy_tb_vars, only: &
       ATMOS_PHY_TB_vars_restart_enddef
    use mod_atmos_phy_cp_vars, only: &
       ATMOS_PHY_CP_vars_restart_enddef
#ifdef _SDM
    use scale_atmos_phy_mp_sdm, only: &
       sd_rest_flg_out, &
       ATMOS_PHY_MP_sdm_restart_enddef
#endif
    implicit none

    !---------------------------------------------------------------------------

#ifdef _SDM
    if( sd_rest_flg_out ) then
       call ATMOS_PHY_MP_sdm_restart_enddef
    endif
#endif

    if ( restart_fid .NE. -1 ) then
       call FILEIO_enddef( restart_fid ) ! [IN]
    endif

    if( ATMOS_sw_dyn )    call ATMOS_DYN_vars_restart_enddef
    if( ATMOS_sw_phy_mp ) call ATMOS_PHY_MP_vars_restart_enddef
    if( ATMOS_sw_phy_ae ) call ATMOS_PHY_AE_vars_restart_enddef
    if( ATMOS_sw_phy_ch ) call ATMOS_PHY_CH_vars_restart_enddef
    if( ATMOS_sw_phy_rd ) call ATMOS_PHY_RD_vars_restart_enddef
    if( ATMOS_sw_phy_sf ) call ATMOS_PHY_SF_vars_restart_enddef
    if( ATMOS_sw_phy_tb ) call ATMOS_PHY_TB_vars_restart_enddef
    if( ATMOS_sw_phy_cp ) call ATMOS_PHY_CP_vars_restart_enddef

    return
  end subroutine ATMOS_vars_restart_enddef

  !-----------------------------------------------------------------------------
  !> Close restart file
  subroutine ATMOS_vars_restart_close
    use scale_fileio, only: &
       FILEIO_close
    use mod_atmos_admin, only: &
       ATMOS_sw_dyn,      &
       ATMOS_sw_phy_mp,   &
       ATMOS_sw_phy_ae,   &
       ATMOS_sw_phy_ch,   &
       ATMOS_sw_phy_rd,   &
       ATMOS_sw_phy_sf,   &
       ATMOS_sw_phy_tb,   &
       ATMOS_sw_phy_cp
    use mod_atmos_dyn_vars, only: &
       ATMOS_DYN_vars_restart_close
    use mod_atmos_phy_mp_vars, only: &
       ATMOS_PHY_MP_vars_restart_close
    use mod_atmos_phy_ae_vars, only: &
       ATMOS_PHY_AE_vars_restart_close
    use mod_atmos_phy_ch_vars, only: &
       ATMOS_PHY_CH_vars_restart_close
    use mod_atmos_phy_rd_vars, only: &
       ATMOS_PHY_RD_vars_restart_close
    use mod_atmos_phy_sf_vars, only: &
       ATMOS_PHY_SF_vars_restart_close
    use mod_atmos_phy_tb_vars, only: &
       ATMOS_PHY_TB_vars_restart_close
    use mod_atmos_phy_cp_vars, only: &
       ATMOS_PHY_CP_vars_restart_close
#ifdef _SDM
    use scale_atmos_phy_mp_sdm, only: &
       sd_rest_flg_out, &
       ATMOS_PHY_MP_sdm_restart_close
#endif
    implicit none

    !---------------------------------------------------------------------------

#ifdef _SDM
    if( sd_rest_flg_out ) then
       call ATMOS_PHY_MP_sdm_restart_close
    endif
#endif

    if ( restart_fid .NE. -1 ) then
       call FILEIO_close( restart_fid ) ! [IN]
       restart_fid = -1

       if ( allocated(VAR_ID) ) deallocate( VAR_ID )
    endif

    if( ATMOS_sw_dyn )    call ATMOS_DYN_vars_restart_close
    if( ATMOS_sw_phy_mp ) call ATMOS_PHY_MP_vars_restart_close
    if( ATMOS_sw_phy_ae ) call ATMOS_PHY_AE_vars_restart_close
    if( ATMOS_sw_phy_ch ) call ATMOS_PHY_CH_vars_restart_close
    if( ATMOS_sw_phy_rd ) call ATMOS_PHY_RD_vars_restart_close
    if( ATMOS_sw_phy_sf ) call ATMOS_PHY_SF_vars_restart_close
    if( ATMOS_sw_phy_tb ) call ATMOS_PHY_TB_vars_restart_close
    if( ATMOS_sw_phy_cp ) call ATMOS_PHY_CP_vars_restart_close

    return
  end subroutine ATMOS_vars_restart_close

  !-----------------------------------------------------------------------------
  !> Define atmospheric variables in restart file
  subroutine ATMOS_vars_restart_def_var
    use scale_fileio, only: &
       FILEIO_def_var
    use mod_atmos_admin, only: &
       ATMOS_sw_dyn,      &
       ATMOS_sw_phy_mp,   &
       ATMOS_sw_phy_ae,   &
       ATMOS_sw_phy_ch,   &
       ATMOS_sw_phy_rd,   &
       ATMOS_sw_phy_sf,   &
       ATMOS_sw_phy_tb,   &
       ATMOS_sw_phy_cp
    use mod_atmos_dyn_vars, only: &
       ATMOS_DYN_vars_restart_def_var
    use mod_atmos_phy_mp_vars, only: &
       ATMOS_PHY_MP_vars_restart_def_var
    use mod_atmos_phy_ae_vars, only: &
       ATMOS_PHY_AE_vars_restart_def_var
    use mod_atmos_phy_ch_vars, only: &
       ATMOS_PHY_CH_vars_restart_def_var
    use mod_atmos_phy_rd_vars, only: &
       ATMOS_PHY_RD_vars_restart_def_var
    use mod_atmos_phy_sf_vars, only: &
       ATMOS_PHY_SF_vars_restart_def_var
    use mod_atmos_phy_tb_vars, only: &
       ATMOS_PHY_TB_vars_restart_def_var
    use mod_atmos_phy_cp_vars, only: &
       ATMOS_PHY_CP_vars_restart_def_var
#ifdef _SDM
    use scale_atmos_phy_mp_sdm, only: &
       sd_rest_flg_out, &
       ATMOS_PHY_MP_sdm_restart_def_var
#endif
    implicit none

    integer iq
    !---------------------------------------------------------------------------

#ifdef _SDM
    if( sd_rest_flg_out ) then
       call ATMOS_PHY_MP_sdm_restart_def_var
    endif
#endif

    if ( restart_fid .NE. -1 ) then

       call FILEIO_def_var( restart_fid, VAR_ID(I_DENS), VAR_NAME(I_DENS), VAR_DESC(I_DENS), & ! [IN]
                            VAR_UNIT(I_DENS), 'ZXY',  ATMOS_RESTART_OUT_DTYPE                ) ! [IN]
       call FILEIO_def_var( restart_fid, VAR_ID(I_MOMZ), VAR_NAME(I_MOMZ), VAR_DESC(I_MOMZ), & ! [IN]
                            VAR_UNIT(I_MOMZ), 'ZHXY', ATMOS_RESTART_OUT_DTYPE                ) ! [IN]
       call FILEIO_def_var( restart_fid, VAR_ID(I_MOMX), VAR_NAME(I_MOMX), VAR_DESC(I_MOMX), & ! [IN]
                            VAR_UNIT(I_MOMX), 'ZXHY', ATMOS_RESTART_OUT_DTYPE                ) ! [IN]
       call FILEIO_def_var( restart_fid, VAR_ID(I_MOMY), VAR_NAME(I_MOMY), VAR_DESC(I_MOMY), & ! [IN]
                            VAR_UNIT(I_MOMY), 'ZXYH', ATMOS_RESTART_OUT_DTYPE                ) ! [IN]
       call FILEIO_def_var( restart_fid, VAR_ID(I_RHOT), VAR_NAME(I_RHOT), VAR_DESC(I_RHOT), & ! [IN]
                            VAR_UNIT(I_RHOT), 'ZXY',  ATMOS_RESTART_OUT_DTYPE                ) ! [IN]

       do iq = 1, QA
          call FILEIO_def_var( restart_fid, VAR_ID(VMAX+iq), TRACER_NAME(iq), TRACER_DESC(iq), & ! [IN]
                               TRACER_UNIT(iq), 'ZXY',  ATMOS_RESTART_OUT_DTYPE       ) ! [IN]
       enddo

    endif

    if( ATMOS_sw_dyn )    call ATMOS_DYN_vars_restart_def_var
    if( ATMOS_sw_phy_mp ) call ATMOS_PHY_MP_vars_restart_def_var
    if( ATMOS_sw_phy_ae ) call ATMOS_PHY_AE_vars_restart_def_var
    if( ATMOS_sw_phy_ch ) call ATMOS_PHY_CH_vars_restart_def_var
    if( ATMOS_sw_phy_rd ) call ATMOS_PHY_RD_vars_restart_def_var
    if( ATMOS_sw_phy_sf ) call ATMOS_PHY_SF_vars_restart_def_var
    if( ATMOS_sw_phy_tb ) call ATMOS_PHY_TB_vars_restart_def_var
    if( ATMOS_sw_phy_cp ) call ATMOS_PHY_CP_vars_restart_def_var

    return
  end subroutine ATMOS_vars_restart_def_var

  !-----------------------------------------------------------------------------
  !> Write restart of atmospheric variables
  subroutine ATMOS_vars_restart_write
    use scale_fileio, only: &
       FILEIO_write_var
    use mod_atmos_admin, only: &
       ATMOS_sw_dyn,      &
       ATMOS_sw_phy_mp,   &
       ATMOS_sw_phy_ae,   &
       ATMOS_sw_phy_ch,   &
       ATMOS_sw_phy_rd,   &
       ATMOS_sw_phy_sf,   &
       ATMOS_sw_phy_tb,   &
       ATMOS_sw_phy_cp
    use mod_atmos_dyn_vars, only: &
       ATMOS_DYN_vars_restart_write
    use mod_atmos_phy_mp_vars, only: &
       ATMOS_PHY_MP_vars_restart_write
    use mod_atmos_phy_ae_vars, only: &
       ATMOS_PHY_AE_vars_restart_write
    use mod_atmos_phy_ch_vars, only: &
       ATMOS_PHY_CH_vars_restart_write
    use mod_atmos_phy_rd_vars, only: &
       ATMOS_PHY_RD_vars_restart_write
    use mod_atmos_phy_sf_vars, only: &
       ATMOS_PHY_SF_vars_restart_write
    use mod_atmos_phy_tb_vars, only: &
       ATMOS_PHY_TB_vars_restart_write
    use mod_atmos_phy_cp_vars, only: &
       ATMOS_PHY_CP_vars_restart_write
#ifdef _SDM
    use scale_atmos_phy_mp_sdm, only: &
       sd_rest_flg_out, &
       ATMOS_PHY_MP_sdm_restart_write
#endif
    implicit none

    integer iq
    !---------------------------------------------------------------------------

#ifdef _SDM
    if( sd_rest_flg_out ) then
       call ATMOS_PHY_MP_sdm_restart_write
    endif
#endif

    if ( restart_fid .NE. -1 ) then

       call ATMOS_vars_fillhalo

       call ATMOS_vars_total

       call FILEIO_write_var( restart_fid, VAR_ID(I_DENS), DENS(:,:,:), VAR_NAME(I_DENS), 'ZXY'  ) ! [IN]
       call FILEIO_write_var( restart_fid, VAR_ID(I_MOMZ), MOMZ(:,:,:), VAR_NAME(I_MOMZ), 'ZHXY' ) ! [IN]
       call FILEIO_write_var( restart_fid, VAR_ID(I_MOMX), MOMX(:,:,:), VAR_NAME(I_MOMX), 'ZXHY' ) ! [IN]
       call FILEIO_write_var( restart_fid, VAR_ID(I_MOMY), MOMY(:,:,:), VAR_NAME(I_MOMY), 'ZXYH' ) ! [IN]
       call FILEIO_write_var( restart_fid, VAR_ID(I_RHOT), RHOT(:,:,:), VAR_NAME(I_RHOT), 'ZXY'  ) ! [IN]

       do iq = 1, QA
          call FILEIO_write_var( restart_fid, VAR_ID(VMAX+iq), QTRC(:,:,:,iq), TRACER_NAME(iq), 'ZXY' ) ! [IN]
       enddo

    endif

    if( ATMOS_sw_dyn )    call ATMOS_DYN_vars_restart_write
    if( ATMOS_sw_phy_mp ) call ATMOS_PHY_MP_vars_restart_write
    if( ATMOS_sw_phy_ae ) call ATMOS_PHY_AE_vars_restart_write
    if( ATMOS_sw_phy_ch ) call ATMOS_PHY_CH_vars_restart_write
    if( ATMOS_sw_phy_rd ) call ATMOS_PHY_RD_vars_restart_write
    if( ATMOS_sw_phy_sf ) call ATMOS_PHY_SF_vars_restart_write
    if( ATMOS_sw_phy_tb ) call ATMOS_PHY_TB_vars_restart_write
    if( ATMOS_sw_phy_cp ) call ATMOS_PHY_CP_vars_restart_write

    return
  end subroutine ATMOS_vars_restart_write

end module mod_atmos_vars
