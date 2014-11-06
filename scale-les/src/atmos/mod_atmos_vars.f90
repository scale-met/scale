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
  public :: ATMOS_vars_history
  public :: ATMOS_vars_total
  public :: ATMOS_vars_diagnostics
  public :: ATMOS_vars_monitor

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical,               public :: ATMOS_RESTART_OUTPUT           = .false.         !< output restart file?
  logical,               public :: ATMOS_RESTART_CHECK            = .false.         !< check value consistency?

  character(len=H_LONG), public :: ATMOS_RESTART_IN_BASENAME      = ''              !< basename of the restart file
  character(len=H_LONG), public :: ATMOS_RESTART_OUT_BASENAME     = ''              !< basename of the output file
  character(len=H_MID),  public :: ATMOS_RESTART_OUT_TITLE        = 'ATMOS restart' !< title    of the output file
  character(len=H_MID),  public :: ATMOS_RESTART_OUT_DTYPE        = 'DEFAULT'       !< REAL4 or REAL8
  logical,               public :: ATMOS_RESTART_IN_ALLOWMISSINGQ = .false.

  character(len=H_LONG), public :: ATMOS_RESTART_CHECK_BASENAME   = 'restart_check'
  real(RP),              public :: ATMOS_RESTART_CHECK_CRITERION  = 1.E-6_RP

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
  logical,                private :: ATMOS_VARS_CHECKRANGE          = .false.

  integer,                private, parameter :: VMAX   = 5       !< number of the variables

  character(len=H_SHORT), private            :: VAR_NAME(VMAX)
  character(len=H_MID),   private            :: VAR_DESC(VMAX)
  character(len=H_SHORT), private            :: VAR_UNIT(VMAX)

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
  integer, private, parameter :: AD_nmax = 51 ! number of diagnostic variables for history output

  integer, private, parameter :: I_W     =  1 ! velocity w at cell center
  integer, private, parameter :: I_U     =  2 ! velocity u at cell center
  integer, private, parameter :: I_V     =  3 ! velocity v at cell center
  integer, private, parameter :: I_POTT  =  4 ! potential temperature

  integer, private, parameter :: I_QDRY  =  5 ! ratio of dry air            to total mass
  integer, private, parameter :: I_QTOT  =  6 ! ratio of total tracer       to total mass
  integer, private, parameter :: I_QHYD  =  7 ! ratio of total hydrometeor  to total mass
  integer, private, parameter :: I_QLIQ  =  8 ! ratio of total liquid water to total mass
  integer, private, parameter :: I_QICE  =  9 ! ratio of total ice    water to total mass

  integer, private, parameter :: I_LWP   = 10 ! liquid water potential temperature
  integer, private, parameter :: I_IWP   = 11 ! relative humidity (liquid+ice)

  integer, private, parameter :: I_RTOT  = 12 ! total gas constant
  integer, private, parameter :: I_CPTOT = 13 ! total heat capacity (constant pressure)
  integer, private, parameter :: I_PRES  = 14 ! pressure
  integer, private, parameter :: I_TEMP  = 15 ! temperature

  integer, private, parameter :: I_POTL  = 16 ! liquid water potential temperature
  integer, private, parameter :: I_RH    = 17 ! relative humidity (liquid+ice)
  integer, private, parameter :: I_RHL   = 18 ! relative humidity against to liquid
  integer, private, parameter :: I_RHI   = 19 ! relative humidity against to ice

  integer, private, parameter :: I_VOR   = 20 ! vertical vorticity
  integer, private, parameter :: I_DIV   = 21 ! divergence
  integer, private, parameter :: I_HDIV  = 22 ! horizontal divergence

  integer, private, parameter :: I_DENS_PRIM = 23 ! prime term of density
  integer, private, parameter :: I_W_PRIM    = 24 ! prime term of w
  integer, private, parameter :: I_U_PRIM    = 25 ! prime term of u
  integer, private, parameter :: I_V_PRIM    = 26 ! prime term of v
  integer, private, parameter :: I_POTT_PRIM = 27 ! prime term of potential temperature
  integer, private, parameter :: I_W_PRIM2   = 28 ! variance of w
  integer, private, parameter :: I_PT_W_PRIM = 29 ! resolved scale heat flux
  integer, private, parameter :: I_W_PRIM3   = 30 ! skewness of w
  integer, private, parameter :: I_TKE_RS    = 31 ! resolved scale TKE

  integer, private, parameter :: I_ENGP  = 32 ! potential energy
  integer, private, parameter :: I_ENGK  = 33 ! kinetic   energy
  integer, private, parameter :: I_ENGI  = 34 ! internal  energy
  integer, private, parameter :: I_ENGT  = 35 ! total     energy

  integer, private, parameter :: I_ENGSFC_SH = 36
  integer, private, parameter :: I_ENGSFC_LH = 37
  integer, private, parameter :: I_ENGSFC_RD = 38
  integer, private, parameter :: I_ENGTOA_RD = 39

  integer, private, parameter :: I_ENGSFC_LW_up = 40
  integer, private, parameter :: I_ENGSFC_LW_dn = 41
  integer, private, parameter :: I_ENGSFC_SW_up = 42
  integer, private, parameter :: I_ENGSFC_SW_dn = 43

  integer, private, parameter :: I_ENGTOA_LW_up = 44
  integer, private, parameter :: I_ENGTOA_LW_dn = 45
  integer, private, parameter :: I_ENGTOA_SW_up = 46
  integer, private, parameter :: I_ENGTOA_SW_dn = 47

  integer, private, parameter :: I_ENGFLXT      = 48

  integer, private, parameter :: I_EVAP         = 49
  integer, private, parameter :: I_PRCP         = 50

  integer, private, parameter :: I_DENS_MEAN    = 51

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
       ATMOS_RESTART_IN_BASENAME,      &
       ATMOS_RESTART_IN_ALLOWMISSINGQ, &
       ATMOS_RESTART_OUTPUT,           &
       ATMOS_RESTART_OUT_BASENAME,     &
       ATMOS_RESTART_OUT_TITLE,        &
       ATMOS_RESTART_OUT_DTYPE,        &
       ATMOS_RESTART_CHECK,            &
       ATMOS_RESTART_CHECK_BASENAME,   &
       ATMOS_RESTART_CHECK_CRITERION,  &
       ATMOS_VARS_CHECKRANGE

    logical :: zinterp ! dummy
    integer :: ierr
    integer :: iv, iq
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[VARS] / Categ[ATMOS] / Origin[SCALE-LES]'

    allocate( DENS(KA,IA,JA)    )
    allocate( MOMZ(KA,IA,JA)    )
    allocate( MOMX(KA,IA,JA)    )
    allocate( MOMY(KA,IA,JA)    )
    allocate( RHOT(KA,IA,JA)    )
    allocate( QTRC(KA,IA,JA,QA) )

    if ( ATMOS_USE_AVERAGE ) then
       allocate( DENS_avw(KA,IA,JA)    )
       allocate( MOMZ_avw(KA,IA,JA)    )
       allocate( MOMX_avw(KA,IA,JA)    )
       allocate( MOMY_avw(KA,IA,JA)    )
       allocate( RHOT_avw(KA,IA,JA)    )
       allocate( QTRC_avw(KA,IA,JA,QA) )

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
    allocate( RHOQ_tp(KA,IA,JA,QA) )

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
                  '*** NO.',5+iq,'|',AQ_NAME(iq),'|', AQ_DESC(iq),'[', AQ_UNIT(iq),']'
    enddo

    if( IO_L ) write(IO_FID_LOG,*)
    if ( ATMOS_RESTART_IN_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart input?  : ', trim(ATMOS_RESTART_IN_BASENAME)
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart input?  : NO'
    endif
    if (       ATMOS_RESTART_OUTPUT             &
         .AND. ATMOS_RESTART_OUT_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart output? : ', trim(ATMOS_RESTART_OUT_BASENAME)
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart output? : NO'
       ATMOS_RESTART_OUTPUT = .false.
    endif

    if( IO_L ) write(IO_FID_LOG,*)
    if (       ATMOS_RESTART_CHECK                &
         .AND. ATMOS_RESTART_CHECK_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Consistency check (for debug)? : YES'
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Consistency check (for debug)? : NO'
       ATMOS_RESTART_CHECK = .false.
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

    allocate( AQ_HIST_id (QA))

    VAR_HIST_id(:) = -1
    AQ_HIST_id (:) = -1
    AD_HIST_id (:) = -1
    AD_MONIT_id(:) = -1
    AD_PREP_sw (:) = -1

    call HIST_reg( VAR_HIST_id(I_DENS), zinterp, VAR_NAME(I_DENS), VAR_DESC(I_DENS), VAR_UNIT(I_DENS), ndim=3 )
    call HIST_reg( VAR_HIST_id(I_MOMZ), zinterp, VAR_NAME(I_MOMZ), VAR_DESC(I_MOMZ), VAR_UNIT(I_MOMZ), ndim=3, zdim='half' )
    call HIST_reg( VAR_HIST_id(I_MOMX), zinterp, VAR_NAME(I_MOMX), VAR_DESC(I_MOMX), VAR_UNIT(I_MOMX), ndim=3, xdim='half' )
    call HIST_reg( VAR_HIST_id(I_MOMY), zinterp, VAR_NAME(I_MOMY), VAR_DESC(I_MOMY), VAR_UNIT(I_MOMY), ndim=3, ydim='half' )
    call HIST_reg( VAR_HIST_id(I_RHOT), zinterp, VAR_NAME(I_RHOT), VAR_DESC(I_RHOT), VAR_UNIT(I_RHOT), ndim=3 )
    do iq = 1, QA
       call HIST_reg( AQ_HIST_id(iq), zinterp, AQ_NAME(iq), AQ_DESC(iq), AQ_UNIT(iq), ndim=3 )
    enddo

    call HIST_reg( AD_HIST_id(I_W   ) , zinterp, 'W',     'velocity w',             'm/s',    ndim=3 )
    call HIST_reg( AD_HIST_id(I_U   ) , zinterp, 'U',     'velocity u',             'm/s',    ndim=3 )
    call HIST_reg( AD_HIST_id(I_V   ) , zinterp, 'V',     'velocity v',             'm/s',    ndim=3 )
    call HIST_reg( AD_HIST_id(I_POTT) , zinterp, 'PT',    'potential temp.',        'K',      ndim=3 )

    call HIST_reg( AD_HIST_id(I_QDRY) , zinterp, 'QDRY',  'dry air',                'kg/kg',  ndim=3 )
    call HIST_reg( AD_HIST_id(I_QTOT) , zinterp, 'QTOT',  'total water',            'kg/kg',  ndim=3 )
    call HIST_reg( AD_HIST_id(I_QHYD) , zinterp, 'QHYD',  'total hydrometeors',     'kg/kg',  ndim=3 )
    call HIST_reg( AD_HIST_id(I_QLIQ) , zinterp, 'QLIQ',  'total liquid water',     'kg/kg',  ndim=3 )
    call HIST_reg( AD_HIST_id(I_QICE) , zinterp, 'QICE',  'total ice water',        'kg/kg',  ndim=3 )

    call HIST_reg( AD_HIST_id(I_LWP)  , zinterp, 'LWP',   'liquid water path',      'g/m2',   ndim=2 )
    call HIST_reg( AD_HIST_id(I_IWP)  , zinterp, 'IWP',   'ice    water path',      'g/m2',   ndim=2 )

    call HIST_reg( AD_HIST_id(I_RTOT) , zinterp, 'RTOT',  'Total gas constant',     'J/kg/K', ndim=3 )
    call HIST_reg( AD_HIST_id(I_CPTOT), zinterp, 'CPTOT', 'Total heat capacity',    'J/kg/K', ndim=3 )
    call HIST_reg( AD_HIST_id(I_PRES) , zinterp, 'PRES',  'pressure',               'Pa',     ndim=3 )
    call HIST_reg( AD_HIST_id(I_TEMP) , zinterp, 'T',     'temperature',            'K',      ndim=3 )

    call HIST_reg( AD_HIST_id(I_POTL) , zinterp, 'LWPT',  'liq. potential temp.',   'K',      ndim=3 )
    call HIST_reg( AD_HIST_id(I_RH)   , zinterp, 'RH',    'relative humidity',      '%',      ndim=3 )
    call HIST_reg( AD_HIST_id(I_RHL)  , zinterp, 'RHL',   'relative humidity(liq)', '%',      ndim=3 )
    call HIST_reg( AD_HIST_id(I_RHI)  , zinterp, 'RHI',   'relative humidity(ice)', '%',      ndim=3 )

    call HIST_reg( AD_HIST_id(I_VOR)  , zinterp, 'VOR',   'vertical vorticity',     '1/s',    ndim=3 )
    call HIST_reg( AD_HIST_id(I_DIV)  , zinterp, 'DIV',   'divergence',             '1/s',    ndim=3 )
    call HIST_reg( AD_HIST_id(I_HDIV) , zinterp, 'HDIV',  'horizontal divergence',  '1/s',    ndim=3 )

    call HIST_reg( AD_HIST_id(I_DENS_PRIM), zinterp, 'DENS_PRIM', 'horiz. deviation of density',    'kg/m3', ndim=3 )
    call HIST_reg( AD_HIST_id(I_W_PRIM   ), zinterp, 'W_PRIM',    'horiz. deviation of w',          'm/s',   ndim=3 )
    call HIST_reg( AD_HIST_id(I_U_PRIM   ), zinterp, 'U_PRIM',    'horiz. deviation of u',          'm/s',   ndim=3 )
    call HIST_reg( AD_HIST_id(I_V_PRIM   ), zinterp, 'V_PRIM',    'horiz. deviation of v',          'm/s',   ndim=3 )
    call HIST_reg( AD_HIST_id(I_POTT_PRIM), zinterp, 'PT_PRIM',   'horiz. deviation of pot. temp.', 'K',     ndim=3 )
    call HIST_reg( AD_HIST_id(I_W_PRIM2  ), zinterp, 'W_PRIM2',   'variance of w',                  'm2/s2', ndim=3 )
    call HIST_reg( AD_HIST_id(I_PT_W_PRIM), zinterp, 'PT_W_PRIM', 'resolved scale heat flux',       'W/s',   ndim=3 )
    call HIST_reg( AD_HIST_id(I_W_PRIM3  ), zinterp, 'W_PRIM3',   'skewness of w',                  'm3/s3', ndim=3 )
    call HIST_reg( AD_HIST_id(I_TKE_RS   ), zinterp, 'TKE_RS',    'resolved scale TKE',             'm2/s2', ndim=3 )

    call HIST_reg( AD_HIST_id(I_ENGT) , zinterp, 'ENGT',  'total energy',           'J/m3',   ndim=3 )
    call HIST_reg( AD_HIST_id(I_ENGP) , zinterp, 'ENGP',  'potential energy',       'J/m3',   ndim=3 )
    call HIST_reg( AD_HIST_id(I_ENGK) , zinterp, 'ENGK',  'kinetic energy',         'J/m3',   ndim=3 )
    call HIST_reg( AD_HIST_id(I_ENGI) , zinterp, 'ENGI',  'internal energy',        'J/m3',   ndim=3 )

    !-----< monitor output setup >-----

    call MONIT_reg( AD_MONIT_id(I_QDRY), 'QDRY', 'dry air mass',     'kg', ndim=3, isflux=.false. )
    call MONIT_reg( AD_MONIT_id(I_QTOT), 'QTOT', 'water mass',       'kg', ndim=3, isflux=.false. )
    call MONIT_reg( AD_MONIT_id(I_EVAP), 'EVAP', 'evaporation',      'kg', ndim=2, isflux=.true.  )
    call MONIT_reg( AD_MONIT_id(I_PRCP), 'PRCP', 'precipitation',    'kg', ndim=2, isflux=.true.  )

    call MONIT_reg( AD_MONIT_id(I_ENGT), 'ENGT', 'total     energy', 'J',  ndim=3, isflux=.false. )
    call MONIT_reg( AD_MONIT_id(I_ENGP), 'ENGP', 'potential energy', 'J',  ndim=3, isflux=.false. )
    call MONIT_reg( AD_MONIT_id(I_ENGK), 'ENGK', 'kinetic   energy', 'J',  ndim=3, isflux=.false. )
    call MONIT_reg( AD_MONIT_id(I_ENGI), 'ENGI', 'internal  energy', 'J',  ndim=3, isflux=.false. )

    call MONIT_reg( AD_MONIT_id(I_ENGFLXT), 'ENGFLXT', 'total energy flux', 'J',  ndim=2, isflux=.true. )

    call MONIT_reg( AD_MONIT_id(I_ENGSFC_SH), 'ENGSFC_SH', 'total     energy', 'J',  ndim=2, isflux=.true. )
    call MONIT_reg( AD_MONIT_id(I_ENGSFC_LH), 'ENGSFC_LH', 'potential energy', 'J',  ndim=2, isflux=.true. )
    call MONIT_reg( AD_MONIT_id(I_ENGSFC_RD), 'ENGSFC_RD', 'kinetic   energy', 'J',  ndim=2, isflux=.true. )
    call MONIT_reg( AD_MONIT_id(I_ENGTOA_RD), 'ENGTOA_RD', 'internal  energy', 'J',  ndim=2, isflux=.true. )

    call MONIT_reg( AD_MONIT_id(I_ENGSFC_LW_up), 'ENGSFC_LW_up', 'total     energy', 'J',  ndim=2, isflux=.true. )
    call MONIT_reg( AD_MONIT_id(I_ENGSFC_LW_dn), 'ENGSFC_LW_dn', 'potential energy', 'J',  ndim=2, isflux=.true. )
    call MONIT_reg( AD_MONIT_id(I_ENGSFC_SW_up), 'ENGSFC_SW_up', 'kinetic   energy', 'J',  ndim=2, isflux=.true. )
    call MONIT_reg( AD_MONIT_id(I_ENGSFC_SW_dn), 'ENGSFC_SW_dn', 'internal  energy', 'J',  ndim=2, isflux=.true. )

    call MONIT_reg( AD_MONIT_id(I_ENGTOA_LW_up), 'ENGTOA_LW_up', 'total     energy', 'J',  ndim=2, isflux=.true. )
    call MONIT_reg( AD_MONIT_id(I_ENGTOA_LW_dn), 'ENGTOA_LW_dn', 'potential energy', 'J',  ndim=2, isflux=.true. )
    call MONIT_reg( AD_MONIT_id(I_ENGTOA_SW_up), 'ENGTOA_SW_up', 'kinetic   energy', 'J',  ndim=2, isflux=.true. )
    call MONIT_reg( AD_MONIT_id(I_ENGTOA_SW_dn), 'ENGTOA_SW_dn', 'internal  energy', 'J',  ndim=2, isflux=.true. )

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
    if (      AD_HIST_id(I_RH)  > 0 &
         .OR. AD_HIST_id(I_RHL) > 0 &
         .OR. AD_HIST_id(I_RHI) > 0 ) then
       AD_PREP_sw(I_QDRY)  = 1
       AD_PREP_sw(I_RTOT)  = 1
       AD_PREP_sw(I_CPTOT) = 1
       AD_PREP_sw(I_PRES)  = 1
       AD_PREP_sw(I_TEMP)  = 1
    endif

    if ( AD_HIST_id (I_VOR) > 0 ) then
       AD_PREP_sw(I_VOR) = 1
    endif

    if ( AD_PREP_sw(I_DIV) > 0 ) then
       AD_PREP_sw(I_HDIV) = 1
    end if

    if ( AD_HIST_id(I_DENS_PRIM) > 0 ) then
       AD_PREP_sw(I_DENS_PRIM) = 1
       AD_PREP_sw(I_DENS_MEAN) = 1
    endif

    if ( AD_HIST_id(I_W_PRIM) > 0 ) then
       AD_PREP_sw(I_W)      = 1
       AD_PREP_sw(I_W_PRIM) = 1
       AD_PREP_sw(I_DENS_MEAN) = 1
    endif

    if ( AD_HIST_id(I_U_PRIM) > 0 ) then
       AD_PREP_sw(I_U)      = 1
       AD_PREP_sw(I_U_PRIM) = 1
       AD_PREP_sw(I_DENS_MEAN) = 1
    endif

    if ( AD_HIST_id(I_V_PRIM) > 0 ) then
       AD_PREP_sw(I_V)      = 1
       AD_PREP_sw(I_V_PRIM) = 1
       AD_PREP_sw(I_DENS_MEAN) = 1
    endif

    if ( AD_HIST_id(I_POTT_PRIM) > 0 ) then
       AD_PREP_sw(I_POTT)      = 1
       AD_PREP_sw(I_POTT_PRIM) = 1
       AD_PREP_sw(I_DENS_MEAN) = 1
    endif

    if ( AD_HIST_id(I_W_PRIM2) > 0 ) then
       AD_PREP_sw(I_W)       = 1
       AD_PREP_sw(I_W_PRIM)  = 1
       AD_PREP_sw(I_W_PRIM2) = 1
       AD_PREP_sw(I_DENS_MEAN) = 1
    endif

    if ( AD_HIST_id(I_PT_W_PRIM) > 0 ) then
       AD_PREP_sw(I_W)         = 1
       AD_PREP_sw(I_W_PRIM)    = 1
       AD_PREP_sw(I_POTT)      = 1
       AD_PREP_sw(I_POTT_PRIM) = 1
       AD_PREP_sw(I_PT_W_PRIM) = 1
       AD_PREP_sw(I_DENS_MEAN) = 1
    endif

    if ( AD_HIST_id(I_W_PRIM3) > 0 ) then
       AD_PREP_sw(I_W)       = 1
       AD_PREP_sw(I_W_PRIM)  = 1
       AD_PREP_sw(I_W_PRIM3) = 1
       AD_PREP_sw(I_DENS_MEAN) = 1
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
    do j  = JS, JE
    do i  = IS, IE
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
    do j  = JS, JE
    do i  = IS, IE
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
  !> Read restart of atmospheric variables
  subroutine ATMOS_vars_restart_read
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       GRAV  => CONST_GRAV
    use scale_fileio, only: &
       FILEIO_read
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

    integer  :: iq
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input restart file (ATMOS) ***'

    if ( ATMOS_RESTART_IN_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(ATMOS_RESTART_IN_BASENAME)

       call FILEIO_read( DENS(:,:,:),                                         & ! [OUT]
                         ATMOS_RESTART_IN_BASENAME, VAR_NAME(1), 'ZXY', step=1 ) ! [IN]
       call FILEIO_read( MOMZ(:,:,:),                                         & ! [OUT]
                         ATMOS_RESTART_IN_BASENAME, VAR_NAME(2), 'ZXY', step=1 ) ! [IN]
       call FILEIO_read( MOMX(:,:,:),                                         & ! [OUT]
                         ATMOS_RESTART_IN_BASENAME, VAR_NAME(3), 'ZXY', step=1 ) ! [IN]
       call FILEIO_read( MOMY(:,:,:),                                         & ! [OUT]
                         ATMOS_RESTART_IN_BASENAME, VAR_NAME(4), 'ZXY', step=1 ) ! [IN]
       call FILEIO_read( RHOT(:,:,:),                                         & ! [OUT]
                         ATMOS_RESTART_IN_BASENAME, VAR_NAME(5), 'ZXY', step=1 ) ! [IN]

       do iq = 1, QA
          call FILEIO_read( QTRC(:,:,:,iq),                                       & ! [OUT]
                            ATMOS_RESTART_IN_BASENAME, AQ_NAME(iq), 'ZXY', step=1 ) ! [IN]
       enddo

       call ATMOS_vars_fillhalo

       call ATMOS_vars_total
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

    if( ATMOS_sw_dyn )    call ATMOS_DYN_vars_restart_read
    if( ATMOS_sw_phy_mp ) call ATMOS_PHY_MP_vars_restart_read
    if( ATMOS_sw_phy_ae ) call ATMOS_PHY_AE_vars_restart_read
    if( ATMOS_sw_phy_ch ) call ATMOS_PHY_CH_vars_restart_read
    if( ATMOS_sw_phy_rd ) call ATMOS_PHY_RD_vars_restart_read
    if( ATMOS_sw_phy_sf ) call ATMOS_PHY_SF_vars_restart_read
    if( ATMOS_sw_phy_tb ) call ATMOS_PHY_TB_vars_restart_read
    if( ATMOS_sw_phy_cp ) call ATMOS_PHY_CP_vars_restart_read

    return
  end subroutine ATMOS_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Write restart of atmospheric variables
  subroutine ATMOS_vars_restart_write
    use scale_time, only: &
       TIME_gettimelabel
    use scale_fileio, only: &
       FILEIO_write
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
    implicit none

    character(len=15)     :: timelabel
    character(len=H_LONG) :: basename

    integer :: iq
    !---------------------------------------------------------------------------

    if ( ATMOS_RESTART_OUT_BASENAME /= '' ) then

       call TIME_gettimelabel( timelabel )
       write(basename,'(A,A,A)') trim(ATMOS_RESTART_OUT_BASENAME), '_', trim(timelabel)

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Output restart file (ATMOS) ***'
       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(basename)

       call ATMOS_vars_fillhalo

       call ATMOS_vars_total

       call FILEIO_write( DENS(:,:,:), basename,                                        ATMOS_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_DENS), VAR_DESC(I_DENS), VAR_UNIT(I_DENS), 'ZXY',  ATMOS_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( MOMZ(:,:,:), basename,                                        ATMOS_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_MOMZ), VAR_DESC(I_MOMZ), VAR_UNIT(I_MOMZ), 'ZHXY', ATMOS_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( MOMX(:,:,:), basename,                                        ATMOS_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_MOMX), VAR_DESC(I_MOMX), VAR_UNIT(I_MOMX), 'ZXHY', ATMOS_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( MOMY(:,:,:), basename,                                        ATMOS_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_MOMY), VAR_DESC(I_MOMY), VAR_UNIT(I_MOMY), 'ZXYH', ATMOS_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( RHOT(:,:,:), basename,                                        ATMOS_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_RHOT), VAR_DESC(I_RHOT), VAR_UNIT(I_RHOT), 'ZXY',  ATMOS_RESTART_OUT_DTYPE  ) ! [IN]

       do iq = 1, QA
          call FILEIO_write( QTRC(:,:,:,iq), basename,                         ATMOS_RESTART_OUT_TITLE, & ! [IN]
                             AQ_NAME(iq), AQ_DESC(iq), AQ_UNIT(iq), 'ZXY',  ATMOS_RESTART_OUT_DTYPE  ) ! [IN]
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

  !-----------------------------------------------------------------------------
  !> Check and compare between last data and sample data
  subroutine ATMOS_vars_restart_check
    use scale_process, only: &
       PRC_myrank
    use gtool_file, only: &
       FileRead
    implicit none

    real(RP) :: DENS_check(KA,IA,JA)    ! Density    [kg/m3]
    real(RP) :: MOMZ_check(KA,IA,JA)    ! momentum z [kg/s/m2]
    real(RP) :: MOMX_check(KA,IA,JA)    ! momentum x [kg/s/m2]
    real(RP) :: MOMY_check(KA,IA,JA)    ! momentum y [kg/s/m2]
    real(RP) :: RHOT_check(KA,IA,JA)    ! DENS * POTT [K*kg/m3]
    real(RP) :: QTRC_check(KA,IA,JA,QA) ! tracer mixing ratio [kg/kg]

    real(RP) :: restart_atmos(KMAX,IMAX,JMAX) !> restart file (no HALO)

    character(len=H_LONG) :: basename

    logical :: datacheck
    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    call PROF_rapstart('Debug')

    write(*,*) 'Compare last Data with ', trim(ATMOS_RESTART_CHECK_BASENAME), 'on PE=', PRC_myrank
    write(*,*) '*** criterion = ', ATMOS_RESTART_CHECK_CRITERION
    datacheck = .true.

    basename = ATMOS_RESTART_CHECK_BASENAME

    call FileRead( restart_atmos(:,:,:), basename, 'DENS', 1, PRC_myrank )
    DENS_check(KS:KE,IS:IE,JS:JE) = restart_atmos(1:KMAX,1:IMAX,1:JMAX)
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

    call FileRead( restart_atmos(:,:,:), basename, 'MOMZ', 1, PRC_myrank )
    MOMZ_check(KS:KE,IS:IE,JS:JE) = restart_atmos(1:KMAX,1:IMAX,1:JMAX)
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

    call FileRead( restart_atmos(:,:,:), basename, 'MOMX', 1, PRC_myrank )
    MOMX_check(KS:KE,IS:IE,JS:JE) = restart_atmos(1:KMAX,1:IMAX,1:JMAX)
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

    call FileRead( restart_atmos(:,:,:), basename, 'MOMY', 1, PRC_myrank )
    MOMY_check(KS:KE,IS:IE,JS:JE) = restart_atmos(1:KMAX,1:IMAX,1:JMAX)
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

    call FileRead( restart_atmos(:,:,:), basename, 'RHOT', 1, PRC_myrank )
    RHOT_check(KS:KE,IS:IE,JS:JE) = restart_atmos(1:KMAX,1:IMAX,1:JMAX)
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
       call FileRead( restart_atmos(:,:,:), basename, AQ_NAME(iq), 1, PRC_myrank )
       QTRC_check(KS:KE,IS:IE,JS:JE,iq) = restart_atmos(1:KMAX,1:IMAX,1:JMAX)
       do k = KS, KE
       do j = JS, JE
       do i = IS, IE
          if ( abs( QTRC(k,i,j,iq)-QTRC_check(k,i,j,iq) ) > ATMOS_RESTART_CHECK_CRITERION ) then
             write(*,*) 'xxx there is the difference  : ', QTRC(k,i,j,iq)-QTRC_check(k,i,j,iq)
             write(*,*) 'xxx at (PE-id,k,i,j,varname) : ', PRC_myrank, k, i, j, AQ_NAME(iq)
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
       LHV   => CONST_LHV,   &
       P00   => CONST_PRE00
    use scale_time, only: &
       TIME_DTSEC
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
       THERMODYN_qd => ATMOS_THERMODYN_qd, &
       CPw => AQ_CP,                       &
       CVw => AQ_CV
    use scale_atmos_saturation, only: &
       SATURATION_dens2qsat_all => ATMOS_SATURATION_dens2qsat_all, &
       SATURATION_dens2qsat_liq => ATMOS_SATURATION_dens2qsat_liq, &
       SATURATION_dens2qsat_ice => ATMOS_SATURATION_dens2qsat_ice
    implicit none

    real(RP) :: QDRY  (KA,IA,JA) ! dry air            [kg/kg]
    real(RP) :: QTOT  (KA,IA,JA) ! total water        [kg/kg]
    real(RP) :: QHYD  (KA,IA,JA) ! total hydrometeor  [kg/kg]
    real(RP) :: QLIQ  (KA,IA,JA) ! total liquid water [kg/kg]
    real(RP) :: QICE  (KA,IA,JA) ! total ice water    [kg/kg]

    real(RP) :: LWP   (IA,JA)    ! liquid water path  [g/m2]
    real(RP) :: IWP   (IA,JA)    ! ice    water path  [g/m2]

    real(RP) :: RTOT  (KA,IA,JA) ! Total gas constant  [J/kg/K]
    real(RP) :: CPTOT (KA,IA,JA) ! Total heat capacity [J/kg/K]
    real(RP) :: CPovCV(KA,IA,JA) ! Cp/Cv

    real(RP) :: POTL  (KA,IA,JA) ! liquid water potential temperature [K]
    real(RP) :: RH    (KA,IA,JA) ! relative humidity (liquid+ice)      [%]
    real(RP) :: RHL   (KA,IA,JA) ! relative humidity against to liquid [%]
    real(RP) :: RHI   (KA,IA,JA) ! relative humidity against to ice    [%]

    real(RP) :: VOR   (KA,IA,JA) ! vertical vorticity    [1/s]
    real(RP) :: DIV   (KA,IA,JA) ! divergence            [1/s]
    real(RP) :: HDIV  (KA,IA,JA) ! horizontal divergence [1/s]

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

    real(RP) :: ENGT  (KA,IA,JA) ! total     energy [J/m3]
    real(RP) :: ENGP  (KA,IA,JA) ! potential energy [J/m3]
    real(RP) :: ENGK  (KA,IA,JA) ! kinetic   energy [J/m3]
    real(RP) :: ENGI  (KA,IA,JA) ! internal  energy [J/m3]

    real(RP) :: QSAT  (KA,IA,JA)
    real(RP) :: UH    (KA,IA,JA)
    real(RP) :: VH    (KA,IA,JA)

    real(RP) :: mean1d(KA)

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
    call HIST_in( DENS(:,:,:), VAR_NAME(I_DENS), VAR_DESC(I_DENS), VAR_UNIT(I_DENS), TIME_DTSEC )
    call HIST_in( MOMZ(:,:,:), VAR_NAME(I_MOMZ), VAR_DESC(I_MOMZ), VAR_UNIT(I_MOMZ), TIME_DTSEC )
    call HIST_in( MOMX(:,:,:), VAR_NAME(I_MOMX), VAR_DESC(I_MOMX), VAR_UNIT(I_MOMX), TIME_DTSEC )
    call HIST_in( MOMY(:,:,:), VAR_NAME(I_MOMY), VAR_DESC(I_MOMY), VAR_UNIT(I_MOMY), TIME_DTSEC )
    call HIST_in( RHOT(:,:,:), VAR_NAME(I_RHOT), VAR_DESC(I_RHOT), VAR_UNIT(I_RHOT), TIME_DTSEC )
    do iq = 1, QA
       call HIST_in( QTRC(:,:,:,iq), AQ_NAME(iq), AQ_DESC(iq), AQ_UNIT(iq), TIME_DTSEC )
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
       call THERMODYN_qd( QDRY(:,:,:),  & ! [OUT]
                          QTRC(:,:,:,:) ) ! [IN]
    endif

    if ( AD_PREP_sw(I_QTOT) > 0 ) then
       do j  = 1, JA
       do i  = 1, IA
       do k  = KS, KE
          QTOT(k,i,j) = 1.0_RP - QDRY(k,i,j)
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_QHYD) > 0 ) then
       QHYD(:,:,:) = 0.0_RP
       do iq = QWS, QWE
         QHYD(:,:,:) = QHYD(:,:,:) + QTRC(:,:,:,iq)
       enddo
       do iq = QIS, QIE
         QHYD(:,:,:) = QHYD(:,:,:) + QTRC(:,:,:,iq)
       enddo
    endif

    if ( AD_PREP_sw(I_QLIQ) > 0 ) then
       QLIQ(:,:,:) = 0.0_RP
       do iq = QWS, QWE
         QLIQ(:,:,:) = QLIQ(:,:,:) + QTRC(:,:,:,iq)
       enddo
    endif

    if ( AD_PREP_sw(I_QICE) > 0 ) then
       QICE(:,:,:) = 0.0_RP
       do iq = QIS, QIE
         QICE(:,:,:) = QICE(:,:,:) + QTRC(:,:,:,iq)
       enddo
    endif

    if ( AD_PREP_sw(I_LWP) > 0 ) then
       do j  = 1, JA
       do i  = 1, IA
          LWP(i,j) = 0.0_RP
          do k  = KS, KE
             LWP(i,j) = LWP(i,j) &
                      + QLIQ(k,i,j) * DENS(k,i,j) * ( REAL_FZ(k,i,j)-REAL_FZ(k-1,i,j) ) * 1.E3_RP ! [kg/m2->g/m2]
          enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_IWP) > 0 ) then
       do j  = 1, JA
       do i  = 1, IA
          IWP(i,j) = 0.0_RP
          do k  = KS, KE
             IWP(i,j) = IWP(i,j) &
                      + QICE(k,i,j) * DENS(k,i,j) * ( REAL_FZ(k,i,j)-REAL_FZ(k-1,i,j) ) * 1.E3_RP ! [kg/m2->g/m2]
          enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_RTOT) > 0 ) then
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j  = 1, JA
       do i  = 1, IA
       do k  = KS, KE
          RTOT (k,i,j) = Rdry * QDRY(k,i,j) + Rvap * QTRC(k,i,j,I_QV)
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_CPTOT) > 0 ) then
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j  = 1, JA
       do i  = 1, IA
       do k  = KS, KE
          CPTOT(k,i,j) = CPdry * QDRY(k,i,j)
       enddo
       enddo
       enddo

       do iq = QQS, QQE
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j  = 1, JA
       do i  = 1, IA
       do k  = KS, KE
          CPTOT(k,i,j) = CPTOT(k,i,j) + QTRC(k,i,j,iq) * CPw(iq)
       enddo
       enddo
       enddo
       enddo

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j  = 1, JA
       do i  = 1, IA
       do k  = KS, KE
          CPovCV(k,i,j) = CPTOT(k,i,j) / ( CPTOT(k,i,j) - RTOT(k,i,j) )
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
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          POTL(k,i,j) = POTT(k,i,j) &
                      - LHV / CPdry * QLIQ(k,i,j) * POTT(k,i,j) / TEMP(k,i,j)
       enddo
       enddo
       enddo
    endif

    if ( AD_HIST_id(I_RH) > 0 ) then
       call SATURATION_dens2qsat_all( QSAT(:,:,:), & ! [OUT]
                                      TEMP(:,:,:), & ! [IN]
                                      DENS(:,:,:)  ) ! [IN]

       do j  = 1, JA
       do i  = 1, IA
       do k  = KS, KE
          RH(k,i,j) = QTRC(k,i,j,I_QV) / QSAT(k,i,j) * 1.0E2_RP
       enddo
       enddo
       enddo
    endif

    if ( AD_HIST_id(I_RHL) > 0 ) then
       call SATURATION_dens2qsat_liq( QSAT(:,:,:), & ! [OUT]
                                      TEMP(:,:,:), & ! [IN]
                                      DENS(:,:,:)  ) ! [IN]

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j  = 1, JA
       do i  = 1, IA
       do k  = KS, KE
          RHL(k,i,j) = QTRC(k,i,j,I_QV) / QSAT(k,i,j) * 1.0E2_RP
       enddo
       enddo
       enddo
    endif

    if ( AD_HIST_id(I_RHI) > 0 ) then
       call SATURATION_dens2qsat_ice( QSAT(:,:,:), & ! [OUT]
                                      TEMP(:,:,:), & ! [IN]
                                      DENS(:,:,:)  ) ! [IN]

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j  = 1, JA
       do i  = 1, IA
       do k  = KS, KE
          RHI(k,i,j) = QTRC(k,i,j,I_QV) / QSAT(k,i,j) * 1.0E2_RP
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_VOR) > 0 ) then
       ! at x, v, layer
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
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
       do j = 2, JA
       do i = 1, IA-1
       do k = KS, KE
          VH(k,i,j) = 0.5_RP * ( MOMY(k,i,j)+MOMY(k,i+1,j)+MOMY(k,i,j-1)+MOMY(k,i+1,j-1) ) &
                             / ( DENS(k,i,j)+DENS(k,i+1,j) )
       enddo
       enddo
       enddo

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
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
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          DIV(k,i,j) = ( MOMZ(k,i,j) - MOMZ(k-1,i  ,j  ) ) * ( REAL_FZ(k,i,j)-REAL_FZ(k-1,i,j) ) &
                     + HDIV(k,i,j)
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_DENS_PRIM) > 0 ) then
       call COMM_horizontal_mean( DENS_MEAN(:), DENS(:,:,:) )
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          DENS_PRIM(k,i,j) = DENS(k,i,j) - DENS_MEAN(k)
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_W_PRIM) > 0 ) then
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          W_PRIM(k,i,j) = W(k,i,j) * DENS(k,i,j)
       enddo
       enddo
       enddo
       call COMM_horizontal_mean( mean1d(:), W_PRIM(:,:,:) )
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          W_PRIM(k,i,j) = ( W_PRIM(k,i,j) - mean1d(k) ) / DENS_MEAN(k)
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_U_PRIM) > 0 ) then
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          U_PRIM(k,i,j) = U(k,i,j) * DENS(k,i,j)
       enddo
       enddo
       enddo
       call COMM_horizontal_mean( mean1d(:), U_PRIM(:,:,:) )
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          U_PRIM(k,i,j) = ( U_PRIM(k,i,j) - mean1d(k) ) / DENS_MEAN(k)
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_V_PRIM) > 0 ) then
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          V_PRIM(k,i,j) = V(k,i,j) * DENS(k,i,j)
       enddo
       enddo
       enddo
       call COMM_horizontal_mean( mean1d(:), V_PRIM(:,:,:) )
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          V_PRIM(k,i,j) = ( V_PRIM(k,i,j) - mean1d(k) ) / DENS_MEAN(k)
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_POTT_PRIM) > 0 ) then
       call COMM_horizontal_mean( mean1d(:), RHOT(:,:,:) )
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          POTT_PRIM(k,i,j) = ( RHOT(k,i,j) - mean1d(k) ) / DENS_MEAN(k)
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_W_PRIM2) > 0 ) then
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          W_PRIM2(k,i,j) = W_PRIM(k,i,j) * W_PRIM(k,i,j)
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_PT_W_PRIM) > 0 ) then
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          PT_W_PRIM(k,i,j) = W_PRIM(k,i,j) * POTT_PRIM(k,i,j) * DENS(k,i,j) * CPdry
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_W_PRIM3) > 0 ) then
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          W_PRIM3(k,i,j) = W_PRIM(k,i,j) * W_PRIM(k,i,j) * W_PRIM(k,i,j)
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_TKE_RS) > 0 ) then
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
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          ENGI(k,i,j) = DENS(k,i,j) * QDRY(k,i,j) * TEMP(k,i,j) * CVdry
       enddo
       enddo
       enddo

       do iq = QQS, QQE
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(3)
       do j  = 1, JA
       do i  = 1, IA
       do k  = KS, KE
          ENGI(k,i,j) = ENGI(k,i,j) &
                      + DENS(k,i,j) * QTRC(k,i,j,iq) * TEMP(k,i,j) * CVw(iq)
       enddo
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_ENGT) > 0 ) then
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j  = 1, JA
       do i  = 1, IA
       do k  = KS, KE
          ENGT(k,i,j) = ENGP(k,i,j) + ENGK(k,i,j) + ENGI(k,i,j)
       enddo
       enddo
       enddo
    endif

    call HIST_in( W    (:,:,:), 'W',     'velocity w',             'm/s',    TIME_DTSEC )
    call HIST_in( U    (:,:,:), 'U',     'velocity u',             'm/s',    TIME_DTSEC )
    call HIST_in( V    (:,:,:), 'V',     'velocity v',             'm/s',    TIME_DTSEC )
    call HIST_in( POTT (:,:,:), 'PT',    'potential temp.',        'K',      TIME_DTSEC )

    call HIST_in( QDRY (:,:,:), 'QDRY',  'dry air',                'kg/kg',  TIME_DTSEC )
    call HIST_in( QTOT (:,:,:), 'QTOT',  'total water',            'kg/kg',  TIME_DTSEC )
    call HIST_in( QHYD (:,:,:), 'QHYD',  'total hydrometeors',     'kg/kg',  TIME_DTSEC )
    call HIST_in( QLIQ (:,:,:), 'QLIQ',  'total liquid water',     'kg/kg',  TIME_DTSEC )
    call HIST_in( QICE (:,:,:), 'QICE',  'total ice water',        'kg/kg',  TIME_DTSEC )

    call HIST_in( LWP  (:,:),   'LWP',   'liquid water path',      'g/m2',   TIME_DTSEC )
    call HIST_in( IWP  (:,:),   'IWP',   'ice    water path',      'g/m2',   TIME_DTSEC )

    call HIST_in( RTOT (:,:,:), 'RTOT',  'Total gas constant',     'J/kg/K', TIME_DTSEC )
    call HIST_in( CPTOT(:,:,:), 'CPTOT', 'Total heat capacity',    'J/kg/K', TIME_DTSEC )
    call HIST_in( PRES (:,:,:), 'PRES',  'pressure',               'Pa',     TIME_DTSEC )
    call HIST_in( TEMP (:,:,:), 'T',     'temperature',            'K',      TIME_DTSEC )

    call HIST_in( POTL (:,:,:), 'LWPT',  'liq. potential temp.',   'K',      TIME_DTSEC )
    call HIST_in( RH   (:,:,:), 'RH',    'relative humidity',      '%',      TIME_DTSEC )
    call HIST_in( RHL  (:,:,:), 'RHL',   'relative humidity(liq)', '%',      TIME_DTSEC )
    call HIST_in( RHI  (:,:,:), 'RHI',   'relative humidity(ice)', '%',      TIME_DTSEC )

    call HIST_in( VOR  (:,:,:), 'VOR',   'vertical vorticity',     '1/s',    TIME_DTSEC )
    call HIST_in( DIV  (:,:,:), 'DIV',   'divergence',             '1/s',    TIME_DTSEC )
    call HIST_in( HDIV (:,:,:), 'HDIV',  'horizontal divergence',  '1/s',    TIME_DTSEC )

    call HIST_in( DENS_PRIM(:,:,:), 'DENS_PRIM', 'horiz. deviation of density',    'kg/m3', TIME_DTSEC )
    call HIST_in( W_PRIM   (:,:,:), 'W_PRIM',    'horiz. deviation of w',          'm/s',   TIME_DTSEC )
    call HIST_in( U_PRIM   (:,:,:), 'U_PRIM',    'horiz. deviation of u',          'm/s',   TIME_DTSEC )
    call HIST_in( V_PRIM   (:,:,:), 'V_PRIM',    'horiz. deviation of v',          'm/s',   TIME_DTSEC )
    call HIST_in( POTT_PRIM(:,:,:), 'PT_PRIM',   'horiz. deviation of pot. temp.', 'K',     TIME_DTSEC )
    call HIST_in( W_PRIM2  (:,:,:), 'W_PRIM2',   'variance of w',                  'm2/s2', TIME_DTSEC )
    call HIST_in( PT_W_PRIM(:,:,:), 'PT_W_PRIM', 'resolved scale heat flux',       'W/s',   TIME_DTSEC )
    call HIST_in( W_PRIM3  (:,:,:), 'W_PRIM3',   'skewness of w',                  'm3/s3', TIME_DTSEC )
    call HIST_in( TKE_RS   (:,:,:), 'TKE_RS',    'resolved scale TKE',             'm2/s2', TIME_DTSEC )

    call HIST_in( ENGT (:,:,:), 'ENGT',  'total energy',           'J/m3',   TIME_DTSEC )
    call HIST_in( ENGP (:,:,:), 'ENGP',  'potential energy',       'J/m3',   TIME_DTSEC )
    call HIST_in( ENGK (:,:,:), 'ENGK',  'kinetic energy',         'J/m3',   TIME_DTSEC )
    call HIST_in( ENGI (:,:,:), 'ENGI',  'internal energy',        'J/m3',   TIME_DTSEC )

    return
  end subroutine ATMOS_vars_history

  !-----------------------------------------------------------------------------
  !> Budget monitor for atmosphere
  subroutine ATMOS_vars_total
    use scale_const, only: &
       GRAV   => CONST_GRAV,   &
       CVdry  => CONST_CVdry
    use scale_grid_real, only: &
       REAL_CZ
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STAT_total
    use scale_atmos_thermodyn, only: &
       THERMODYN_qd        => ATMOS_THERMODYN_qd,        &
       THERMODYN_temp_pres => ATMOS_THERMODYN_temp_pres, &
       CVw => AQ_CV
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

          call STAT_total( total, RHOQ(:,:,:), AQ_NAME(iq) )
       enddo

       call THERMODYN_qd( QDRY(:,:,:),  & ! [OUT]
                          QTRC(:,:,:,:) ) ! [IN]

       call THERMODYN_temp_pres( TEMP(:,:,:),  & ! [OUT]
                                 PRES(:,:,:),  & ! [OUT]
                                 DENS(:,:,:),  & ! [IN]
                                 RHOT(:,:,:),  & ! [IN]
                                 QTRC(:,:,:,:) ) ! [IN]

       RHOQ(:,:,:) = DENS(:,:,:) * QDRY (:,:,:)

       call STAT_total( total, RHOQ(:,:,:), 'QDRY' )

       RHOQ(:,:,:) = DENS(:,:,:) * ( 1.0_RP - QDRY (:,:,:) ) ! Qtotal

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
          do iq = QQS, QQE
             ENGI(k,i,j) = ENGI(k,i,j) &
                         + DENS(k,i,j) * QTRC(k,i,j,iq) * TEMP(k,i,j) * CVw(iq)
          enddo

          ENGT(k,i,j) = ENGP(k,i,j) + ENGK(k,i,j) + ENGI(k,i,j)
       enddo
       enddo
       enddo

       call STAT_total( total, ENGT(:,:,:), 'ENGT' )
       call STAT_total( total, ENGP(:,:,:), 'ENGP' )
       call STAT_total( total, ENGK(:,:,:), 'ENGK' )
       call STAT_total( total, ENGI(:,:,:), 'ENGI' )

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

    call THERMODYN_temp_pres( TEMP(:,:,:),  & ! [OUT]
                              PRES(:,:,:),  & ! [OUT]
                              DENS(:,:,:),  & ! [IN]
                              RHOT(:,:,:),  & ! [IN]
                              QTRC(:,:,:,:) ) ! [IN]

    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       W(k,i,j) = 0.5_RP * ( MOMZ(k-1,i,j)+MOMZ(k,i,j) ) / DENS(k,i,j)
    enddo
    enddo
    enddo

    do j = 1, JA
    do i = 2, IA
    do k = KS, KE
       U(k,i,j) = 0.5_RP * ( MOMX(k,i-1,j)+MOMX(k,i,j) ) / DENS(k,i,j)
    enddo
    enddo
    enddo
    do j = 1, JA
    do k = KS, KE
       U(k,1,j) = MOMX(k,1,j) / DENS(k,1,j)
    enddo
    enddo

    do j = 2, JA
    do i = 1, IA
    do k = KS, KE
       V(k,i,j) = 0.5_RP * ( MOMY(k,i,j-1)+MOMY(k,i,j) ) / DENS(k,i,j)
    enddo
    enddo
    enddo
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
    use scale_const, only: &
       GRAV   => CONST_GRAV,   &
       CVdry  => CONST_CVdry
    use scale_grid_real, only: &
       REAL_CZ
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STAT_total,            &
       STAT_detail
    use scale_monitor, only: &
       MONIT_put, &
       MONIT_in
    use scale_atmos_thermodyn, only: &
       THERMODYN_qd        => ATMOS_THERMODYN_qd,        &
       THERMODYN_temp_pres => ATMOS_THERMODYN_temp_pres, &
       CVw => AQ_CV
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
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          RHOQ(k,i,j) = DENS(k,i,j) * QTRC(k,i,j,iq)
       enddo
       enddo
       enddo

       call MONIT_in( RHOQ(:,:,:), AQ_NAME(iq), AQ_DESC(iq), AQ_UNIT(iq), ndim=3, isflux=.false. )
    enddo

    ! total dry airmass

    call THERMODYN_qd( QDRY(:,:,:),  & ! [OUT]
                       QTRC(:,:,:,:) ) ! [IN]

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
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
    do j = JS, JE
    do i = IS, IE
       PRCP(i,j) = SFLX_rain(i,j) + SFLX_snow(i,j)
    enddo
    enddo
    call MONIT_put( AD_MONIT_id(I_PRCP), PRCP(:,:) )

    !##### Energy Budget #####

    call THERMODYN_temp_pres( TEMP(:,:,:),  & ! [OUT]
                              PRES(:,:,:),  & ! [OUT]
                              DENS(:,:,:),  & ! [IN]
                              RHOT(:,:,:),  & ! [IN]
                              QTRC(:,:,:,:) ) ! [IN]

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       ENGP(k,i,j) = DENS(k,i,j) * GRAV * REAL_CZ(k,i,j)

       ENGK(k,i,j) = 0.5_RP * DENS(k,i,j) * ( W(k,i,j)**2 &
                                            + U(k,i,j)**2 &
                                            + V(k,i,j)**2 )

       ENGI(k,i,j) = DENS(k,i,j) * QDRY(k,i,j) * TEMP(k,i,j) * CVdry
       do iq = QQS, QQE
          ENGI(k,i,j) = ENGI(k,i,j) &
                      + DENS(k,i,j) * QTRC(k,i,j,iq) * TEMP(k,i,j) * CVw(iq)
       enddo
    enddo
    enddo
    enddo

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       ENGT(k,i,j) = ENGP(k,i,j) + ENGK(k,i,j) + ENGI(k,i,j)
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       SFLX_RD_net(i,j) = ( SFLX_LW_up(i,j) - SFLX_LW_dn(i,j) ) &
                        + ( SFLX_SW_up(i,j) - SFLX_SW_dn(i,j) )

       TOAFLX_RD_net(i,j) = ( TOAFLX_LW_up(i,j) - TOAFLX_LW_dn(i,j) ) &
                          + ( TOAFLX_SW_up(i,j) - TOAFLX_SW_dn(i,j) )
    enddo
    enddo

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
       ENGFLXT(i,j) = SFLX_SH(i,j) + SFLX_LH(i,j) &
                    + SFLX_RD_net(i,j) - TOAFLX_RD_net(i,j)
    enddo
    enddo

    call MONIT_put( AD_MONIT_id(I_ENGT), ENGT(:,:,:) )
    call MONIT_put( AD_MONIT_id(I_ENGP), ENGP(:,:,:) )
    call MONIT_put( AD_MONIT_id(I_ENGK), ENGK(:,:,:) )
    call MONIT_put( AD_MONIT_id(I_ENGI), ENGI(:,:,:) )

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
       WORK(:,:,:,1) = W(:,:,:)
       WORK(:,:,:,2) = U(:,:,:)
       WORK(:,:,:,3) = V(:,:,:)

       WNAME(1) = "W"
       WNAME(2) = "U"
       WNAME(3) = "V"

       call STAT_detail( WORK(:,:,:,:), WNAME(:) )
    endif

    return
  end subroutine ATMOS_vars_monitor

end module mod_atmos_vars
