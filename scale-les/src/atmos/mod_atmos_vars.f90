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

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
# include "scalelib.h"

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  character(len=H_SHORT), public, save :: ATMOS_DYN_TYPE    = 'NONE'
  character(len=H_SHORT), public, save :: ATMOS_PHY_SF_TYPE = 'NONE'
  character(len=H_SHORT), public, save :: ATMOS_PHY_TB_TYPE = 'NONE'
  character(len=H_SHORT), public, save :: ATMOS_PHY_MP_TYPE = 'NONE'
  character(len=H_SHORT), public, save :: ATMOS_PHY_RD_TYPE = 'NONE'
  character(len=H_SHORT), public, save :: ATMOS_PHY_AE_TYPE = 'NONE'
  character(len=H_SHORT), public, save :: ATMOS_PHY_CH_TYPE = 'NONE'

  logical,                public, save :: ATMOS_sw_dyn
  logical,                public, save :: ATMOS_sw_phy_sf
  logical,                public, save :: ATMOS_sw_phy_tb
  logical,                public, save :: ATMOS_sw_phy_mp
  logical,                public, save :: ATMOS_sw_phy_rd
  logical,                public, save :: ATMOS_sw_phy_ae
  logical,                public, save :: ATMOS_sw_restart
  logical,                public, save :: ATMOS_sw_check

  logical,                public, save :: ATMOS_USE_AVERAGE = .false.

  ! prognostic variables
  real(RP), public, target, allocatable :: DENS(:,:,:)   ! Density     [kg/m3]
  real(RP), public, target, allocatable :: MOMZ(:,:,:)   ! momentum z  [kg/s/m2]
  real(RP), public, target, allocatable :: MOMX(:,:,:)   ! momentum x  [kg/s/m2]
  real(RP), public, target, allocatable :: MOMY(:,:,:)   ! momentum y  [kg/s/m2]
  real(RP), public, target, allocatable :: RHOT(:,:,:)   ! DENS * POTT [K*kg/m3]
  real(RP), public, target, allocatable :: QTRC(:,:,:,:) ! ratio of mass of tracer to total mass[kg/kg]

  real(RP), public, target, allocatable :: DENS_avw(:,:,:)
  real(RP), public, target, allocatable :: MOMZ_avw(:,:,:)
  real(RP), public, target, allocatable :: MOMX_avw(:,:,:)
  real(RP), public, target, allocatable :: MOMY_avw(:,:,:)
  real(RP), public, target, allocatable :: RHOT_avw(:,:,:)
  real(RP), public, target, allocatable :: QTRC_avw(:,:,:,:)

  real(RP), public, pointer :: DENS_av(:,:,:)
  real(RP), public, pointer :: MOMZ_av(:,:,:)
  real(RP), public, pointer :: MOMX_av(:,:,:)
  real(RP), public, pointer :: MOMY_av(:,:,:)
  real(RP), public, pointer :: RHOT_av(:,:,:)
  real(RP), public, pointer :: QTRC_av(:,:,:,:)

  ! tendency by physical processes
  real(RP), public, allocatable :: DENS_tp(:,:,:)
  real(RP), public, allocatable :: MOMZ_tp(:,:,:)
  real(RP), public, allocatable :: MOMX_tp(:,:,:)
  real(RP), public, allocatable :: MOMY_tp(:,:,:)
  real(RP), public, allocatable :: RHOT_tp(:,:,:)
  real(RP), public, allocatable :: QTRC_tp(:,:,:,:)

  character(len=H_SHORT), public, save :: AP_NAME(5)
  character(len=H_MID),   public, save :: AP_DESC(5)
  character(len=H_SHORT), public, save :: AP_UNIT(5)

  data AP_NAME / 'DENS', &
                 'MOMZ', &
                 'MOMX', &
                 'MOMY', &
                 'RHOT'  /
  data AP_DESC / 'density',    &
                 'momentum z', &
                 'momentum x', &
                 'momentum y', &
                 'rho * theta' /
  data AP_UNIT / 'kg/m3',   &
                 'kg/m2/s', &
                 'kg/m2/s', &
                 'kg/m2/s', &
                 'kg/m3*K'  /

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical,               private, save :: ATMOS_RESTART_OUTPUT           = .false.
  character(len=H_LONG), private, save :: ATMOS_RESTART_IN_BASENAME      = 'restart_in'
  character(len=H_LONG), private, save :: ATMOS_RESTART_OUT_BASENAME     = 'restart_out'
  character(len=H_MID),  private, save :: ATMOS_RESTART_OUT_TITLE        = 'SCALE-LES PROGNOSTIC VARS.'
  logical,               private, save :: ATMOS_RESTART_IN_ALLOWMISSINGQ = .false.

  logical,               private, save :: ATMOS_RESTART_CHECK            = .false.
  character(len=H_LONG), private, save :: ATMOS_RESTART_CHECK_BASENAME   = 'restart_check'
  real(RP),              private, save :: ATMOS_RESTART_CHECK_CRITERION  = 1.E-6_RP

  logical,               private, save :: ATMOS_VARS_CHECKRANGE          = .false.

  ! history output of prognostic variables
  integer, private              :: AP_HIST_id(5)
  integer, private, allocatable :: AQ_HIST_id(:)

  ! history & monitor output of diagnostic variables
  integer, private, parameter :: AD_nmax = 26 ! number of diagnostic variables for history output

  integer, private, parameter :: I_VELZ  =  1 ! velocity w at cell center
  integer, private, parameter :: I_VELX  =  2 ! velocity u at cell center
  integer, private, parameter :: I_VELY  =  3 ! velocity v at cell center
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

  integer, private, parameter :: I_ENGP  = 23 ! potential energy
  integer, private, parameter :: I_ENGK  = 24 ! kinetic   energy
  integer, private, parameter :: I_ENGI  = 25 ! internal  energy
  integer, private, parameter :: I_ENGT  = 26 ! total     energy

  integer, private, save      :: AD_HIST_id (AD_nmax)
  integer, private, save      :: AD_PREP_sw (AD_nmax)
  integer, private, save      :: AD_MONIT_id(AD_nmax)

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
    implicit none

    NAMELIST / PARAM_ATMOS / &
       ATMOS_DYN_TYPE, &
       ATMOS_PHY_SF_TYPE, &
       ATMOS_PHY_TB_TYPE, &
       ATMOS_PHY_MP_TYPE, &
       ATMOS_PHY_RD_TYPE, &
       ATMOS_PHY_AE_TYPE, &
       ATMOS_PHY_CH_TYPE

    NAMELIST / PARAM_ATMOS_VARS / &
       ATMOS_RESTART_IN_BASENAME,      &
       ATMOS_RESTART_IN_ALLOWMISSINGQ, &
       ATMOS_RESTART_OUTPUT,           &
       ATMOS_RESTART_OUT_BASENAME,     &
       ATMOS_RESTART_CHECK,            &
       ATMOS_RESTART_CHECK_BASENAME,   &
       ATMOS_RESTART_CHECK_CRITERION,  &
       ATMOS_VARS_CHECKRANGE

    logical :: zinterp ! dummy
    integer :: ierr
    integer :: ip, iq
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Variables]/Categ[ATMOS]'

    allocate( AQ_HIST_id(QA) )

    allocate( DENS(KA,IA,JA) )
    allocate( MOMZ(KA,IA,JA) )
    allocate( MOMX(KA,IA,JA) )
    allocate( MOMY(KA,IA,JA) )
    allocate( RHOT(KA,IA,JA) )
    allocate( QTRC(KA,IA,JA,QA) )

    allocate( DENS_tp(KA,IA,JA) )
    allocate( MOMZ_tp(KA,IA,JA) )
    allocate( MOMX_tp(KA,IA,JA) )
    allocate( MOMY_tp(KA,IA,JA) )
    allocate( RHOT_tp(KA,IA,JA) )
    allocate( QTRC_tp(KA,IA,JA,QA) )


    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS)

    !-----< module component check >-----

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** [ATMOS] selected components'

    if( IO_L ) write(IO_FID_LOG,*) 'Dynamics...'

    if ( ATMOS_DYN_TYPE /= 'OFF' .AND. ATMOS_DYN_TYPE /= 'NONE' ) then
       if( IO_L ) write(IO_FID_LOG,*) '  Dynamical core   : ON'
       if( IO_L ) write(IO_FID_LOG,*) '  Tracer advection : ON'
       ATMOS_sw_dyn = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '  Dynamical core   : OFF'
       if( IO_L ) write(IO_FID_LOG,*) '  Tracer advection : OFF'
       ATMOS_sw_dyn = .false.
    endif

    if( IO_L ) write(IO_FID_LOG,*) 'Physics...'

    if ( ATMOS_PHY_SF_TYPE /= 'OFF' .AND. ATMOS_PHY_SF_TYPE /= 'NONE' ) then
       if( IO_L ) write(IO_FID_LOG,*) '  Surface Flux : ON'
       ATMOS_sw_phy_sf = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '  Surface Flux : OFF'
       ATMOS_sw_phy_sf = .false.
    endif
    if ( ATMOS_PHY_TB_TYPE /= 'OFF' .AND. ATMOS_PHY_TB_TYPE /= 'NONE' ) then
       if( IO_L ) write(IO_FID_LOG,*) '  Sub-grid Turbulence : ON'
       ATMOS_sw_phy_tb = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '  Sub-grid Turbulence : OFF'
       ATMOS_sw_phy_tb = .false.
    endif
    if ( ATMOS_PHY_MP_TYPE /= 'OFF' .AND. ATMOS_PHY_MP_TYPE /= 'NONE' ) then
       if( IO_L ) write(IO_FID_LOG,*) '  Cloud Microphysics  : ON'
       ATMOS_sw_phy_mp = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '  Cloud Microphysics  : OFF'
       ATMOS_sw_phy_mp = .false.
    endif
    if ( ATMOS_PHY_RD_TYPE /= 'OFF' .AND. ATMOS_PHY_RD_TYPE /= 'NONE' ) then
       if( IO_L ) write(IO_FID_LOG,*) '  Radiative transfer  : ON'
       ATMOS_sw_phy_rd = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '  Radiative transfer  : OFF'
       ATMOS_sw_phy_rd = .false.
    endif
    if ( ATMOS_PHY_AE_TYPE /= 'OFF' .AND. ATMOS_PHY_AE_TYPE /= 'NONE' ) then
       if( IO_L ) write(IO_FID_LOG,*) '  Aerosol  : ON'
       ATMOS_sw_phy_ae = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '  Aerosol  : OFF'
       ATMOS_sw_phy_ae = .false.
    endif

    !-----< prognostic variable list check >-----

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[ATMOS VARS]/Categ[ATMOS]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_VARS,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_VARS. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_VARS)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** [ATMOS] prognostic variables'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A8,A,A32,3(A))') &
               '***       |',' VARNAME','|', 'DESCRIPTION                     ','[', 'UNIT            ',']'
    do ip = 1, 5
       if( IO_L ) write(IO_FID_LOG,'(1x,A,i3,A,A8,A,A32,3(A))') &
                  '*** NO.',ip,'|',trim(AP_NAME(ip)),'|', AP_DESC(ip),'[', AP_UNIT(ip),']'
    enddo
    do iq = 1, QA
       if( IO_L ) write(IO_FID_LOG,'(1x,A,i3,A,A8,A,A32,3(A))')  &
                  '*** NO.',5+iq,'|',trim(AQ_NAME(iq)),'|', AQ_DESC(iq),'[', AQ_UNIT(iq),']'
    enddo

    !-----< restart output & consistency check >-----

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) 'Output...'
    if ( ATMOS_RESTART_OUTPUT ) then
       if( IO_L ) write(IO_FID_LOG,*) '  Atmos restart output : YES'
       ATMOS_sw_restart = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '  Atmos restart output : NO'
       ATMOS_sw_restart = .false.
    endif
    if ( ATMOS_RESTART_CHECK ) then
       if( IO_L ) write(IO_FID_LOG,*) '  Data check : YES'
       ATMOS_sw_check = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '  Data check : NO'
       ATMOS_sw_check = .false.
    endif
    if( IO_L ) write(IO_FID_LOG,*)

    !-----< use average? >-----

    if ( ATMOS_USE_AVERAGE ) then
       if( IO_L ) write(IO_FID_LOG,*) '  Atmos use average : YES'
       allocate( DENS_avw(KA,IA,JA) )
       allocate( MOMZ_avw(KA,IA,JA) )
       allocate( MOMX_avw(KA,IA,JA) )
       allocate( MOMY_avw(KA,IA,JA) )
       allocate( RHOT_avw(KA,IA,JA) )
       allocate( QTRC_avw(KA,IA,JA,QA) )

       DENS_av => DENS_avw
       MOMZ_av => MOMZ_avw
       MOMX_av => MOMX_avw
       MOMY_av => MOMY_avw
       RHOT_av => RHOT_avw
       QTRC_av => QTRC_avw
    else
       if( IO_L ) write(IO_FID_LOG,*) '  Atmos use average : NO'
       DENS_av => DENS
       MOMZ_av => MOMZ
       MOMX_av => MOMX
       MOMY_av => MOMY
       RHOT_av => RHOT
       QTRC_av => QTRC
    end if

    !-----< history output setup: general set >-----

    AP_HIST_id    (:) = -1
    AQ_HIST_id    (:) = -1
    AD_HIST_id (:) = -1
    AD_MONIT_id(:) = -1
    AD_PREP_sw (:) = -1

    call HIST_reg( AP_HIST_id(I_DENS), zinterp, AP_NAME(I_DENS), AP_DESC(I_DENS), AP_UNIT(I_DENS), ndim=3 )
    call HIST_reg( AP_HIST_id(I_MOMZ), zinterp, AP_NAME(I_MOMZ), AP_DESC(I_MOMZ), AP_UNIT(I_MOMZ), ndim=3, zdim='half' )
    call HIST_reg( AP_HIST_id(I_MOMX), zinterp, AP_NAME(I_MOMX), AP_DESC(I_MOMX), AP_UNIT(I_MOMX), ndim=3, xdim='half' )
    call HIST_reg( AP_HIST_id(I_MOMY), zinterp, AP_NAME(I_MOMY), AP_DESC(I_MOMY), AP_UNIT(I_MOMY), ndim=3, ydim='half' )
    call HIST_reg( AP_HIST_id(I_RHOT), zinterp, AP_NAME(I_RHOT), AP_DESC(I_RHOT), AP_UNIT(I_RHOT), ndim=3 )
    do iq = 1, QA
       call HIST_reg( AQ_HIST_id(iq), zinterp, AQ_NAME(iq), AQ_DESC(iq), AQ_UNIT(iq), ndim=3 )
    enddo

    call HIST_reg( AD_HIST_id(I_VELZ) , zinterp, 'W',     'velocity w',             'm/s',    ndim=3 )
    call HIST_reg( AD_HIST_id(I_VELX) , zinterp, 'U',     'velocity u',             'm/s',    ndim=3 )
    call HIST_reg( AD_HIST_id(I_VELY) , zinterp, 'V',     'velocity v',             'm/s',    ndim=3 )
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

    call HIST_reg( AD_HIST_id(I_ENGT) , zinterp, 'ENGT',  'total energy',           'J/m3',   ndim=3 )
    call HIST_reg( AD_HIST_id(I_ENGP) , zinterp, 'ENGP',  'potential energy',       'J/m3',   ndim=3 )
    call HIST_reg( AD_HIST_id(I_ENGK) , zinterp, 'ENGK',  'kinetic energy',         'J/m3',   ndim=3 )
    call HIST_reg( AD_HIST_id(I_ENGI) , zinterp, 'ENGI',  'internal energy',        'J/m3',   ndim=3 )

    !-----< monitor output setup >-----

    call MONIT_reg( AD_MONIT_id(I_QDRY), 'QDRY', 'dry air mass',     'kg', ndim=3 )
    call MONIT_reg( AD_MONIT_id(I_QTOT), 'QTOT', 'water mass',       'kg', ndim=3 )
    call MONIT_reg( AD_MONIT_id(I_ENGT), 'ENGT', 'total     energy', 'J',  ndim=3 )
    call MONIT_reg( AD_MONIT_id(I_ENGP), 'ENGP', 'potential energy', 'J',  ndim=3 )
    call MONIT_reg( AD_MONIT_id(I_ENGK), 'ENGK', 'kinetic   energy', 'J',  ndim=3 )
    call MONIT_reg( AD_MONIT_id(I_ENGI), 'ENGI', 'internal  energy', 'J',  ndim=3 )

    if ( AD_HIST_id(I_VELZ) > 0 ) then
       AD_PREP_sw(I_VELZ) = 1
    endif
    if ( AD_HIST_id(I_VELX) > 0 ) then
       AD_PREP_sw(I_VELX) = 1
    endif
    if ( AD_HIST_id(I_VELY) > 0 ) then
       AD_PREP_sw(I_VELY) = 1
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

    if (      AD_HIST_id (I_ENGP) > 0 &
         .OR. AD_MONIT_id(I_ENGP) > 0 ) then
       AD_PREP_sw(I_ENGP) = 1
    endif
    if (      AD_HIST_id (I_ENGK) > 0 &
         .OR. AD_MONIT_id(I_ENGK) > 0 ) then
       AD_PREP_sw(I_VELZ) = 1
       AD_PREP_sw(I_VELX) = 1
       AD_PREP_sw(I_VELY) = 1
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
       AD_PREP_sw(I_VELZ)  = 1
       AD_PREP_sw(I_VELX)  = 1
       AD_PREP_sw(I_VELY)  = 1
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
  !> fill HALO region of atmospheric variables
  subroutine ATMOS_vars_fillhalo
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    integer :: i, j, iq
    !---------------------------------------------------------------------------

    ! fill KHALO
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

    ! fill IHALO & JHALO
    call COMM_vars8( DENS(:,:,:), 1 )
    call COMM_vars8( MOMZ(:,:,:), 2 )
    call COMM_vars8( MOMX(:,:,:), 3 )
    call COMM_vars8( MOMY(:,:,:), 4 )
    call COMM_vars8( RHOT(:,:,:), 5 )
    call COMM_wait ( DENS(:,:,:), 1 )
    call COMM_wait ( MOMZ(:,:,:), 2 )
    call COMM_wait ( MOMX(:,:,:), 3 )
    call COMM_wait ( MOMY(:,:,:), 4 )
    call COMM_wait ( RHOT(:,:,:), 5 )

    do iq = 1, QA
       call COMM_vars8( QTRC(:,:,:,iq), iq )
    enddo
    do iq = 1, QA
       call COMM_wait ( QTRC(:,:,:,iq), iq )
    enddo

    return
  end subroutine ATMOS_vars_fillhalo

  !-----------------------------------------------------------------------------
  !> Read restart of atmospheric variables
  subroutine ATMOS_vars_restart_read
    use gtool_file, only: &
       FileRead
    use scale_process, only: &
       PRC_myrank
    use scale_const, only: &
       GRAV  => CONST_GRAV,  &
       CVdry => CONST_CVdry
    use scale_grid, only: &
       CZ => GRID_CZ
    use scale_monitor, only: &
       MONIT_put,  &
       MONIT_in,   &
       MONIT_write
    use scale_atmos_thermodyn, only: &
       THERMODYN_qd        => ATMOS_THERMODYN_qd,        &
       THERMODYN_temp_pres => ATMOS_THERMODYN_temp_pres, &
       CVw => AQ_CV
    implicit none

    real(RP) :: restart_atmos(KMAX,IMAX,JMAX) !> restart file (no HALO)

    character(len=H_LONG) :: basename

    real(RP) :: VELZ  (KA,IA,JA) ! velocity w at cell center [m/s]
    real(RP) :: VELX  (KA,IA,JA) ! velocity u at cell center [m/s]
    real(RP) :: VELY  (KA,IA,JA) ! velocity v at cell center [m/s]

    real(RP) :: QDRY  (KA,IA,JA) ! dry air     [kg/kg]
    real(RP) :: PRES  (KA,IA,JA) ! pressure    [Pa]
    real(RP) :: TEMP  (KA,IA,JA) ! temperature [K]

    real(RP) :: ENGT  (KA,IA,JA) ! total     energy [J/m3]
    real(RP) :: ENGP  (KA,IA,JA) ! potential energy [J/m3]
    real(RP) :: ENGK  (KA,IA,JA) ! kinetic   energy [J/m3]
    real(RP) :: ENGI  (KA,IA,JA) ! internal  energy [J/m3]

    real(RP) :: RHOQ  (KA,IA,JA)

    integer :: i, j, k, iq
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input restart file (atmos) ***'

    call PROF_rapstart('FILE I NetCDF')

    basename = ATMOS_RESTART_IN_BASENAME

    call FileRead( restart_atmos(:,:,:), basename, 'DENS', 1, PRC_myrank )
    DENS(KS:KE,IS:IE,JS:JE) = restart_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FileRead( restart_atmos(:,:,:), basename, 'MOMZ', 1, PRC_myrank )
    MOMZ(KS:KE,IS:IE,JS:JE) = restart_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FileRead( restart_atmos(:,:,:), basename, 'MOMX', 1, PRC_myrank )
    MOMX(KS:KE,IS:IE,JS:JE) = restart_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FileRead( restart_atmos(:,:,:), basename, 'MOMY', 1, PRC_myrank )
    MOMY(KS:KE,IS:IE,JS:JE) = restart_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FileRead( restart_atmos(:,:,:), basename, 'RHOT', 1, PRC_myrank )
    RHOT(KS:KE,IS:IE,JS:JE) = restart_atmos(1:KMAX,1:IMAX,1:JMAX)

    do iq = 1, QA
       call FileRead( restart_atmos(:,:,:), basename, AQ_NAME(iq), 1, PRC_myrank )
       QTRC(KS:KE,IS:IE,JS:JE,iq) = restart_atmos(1:KMAX,1:IMAX,1:JMAX)
    enddo

    call PROF_rapend  ('FILE I NetCDF')

    ! fill halo
    call ATMOS_vars_fillhalo

    ! check total (optional)
    call ATMOS_vars_total

    if ( ATMOS_USE_AVERAGE ) then
       DENS_av(:,:,:) = DENS(:,:,:)
       MOMZ_av(:,:,:) = MOMZ(:,:,:)
       MOMX_av(:,:,:) = MOMX(:,:,:)
       MOMY_av(:,:,:) = MOMY(:,:,:)
       RHOT_av(:,:,:) = RHOT(:,:,:)

       QTRC_av(:,:,:,:) = QTRC(:,:,:,:)
    endif

    ! first monitor output
    call MONIT_in( DENS(:,:,:), AP_NAME(I_DENS), AP_DESC(I_DENS), AP_UNIT(I_DENS), ndim=3 )
    call MONIT_in( MOMZ(:,:,:), AP_NAME(I_MOMZ), AP_DESC(I_MOMZ), AP_UNIT(I_MOMZ), ndim=3 )
    call MONIT_in( MOMX(:,:,:), AP_NAME(I_MOMX), AP_DESC(I_MOMX), AP_UNIT(I_MOMX), ndim=3 )
    call MONIT_in( MOMY(:,:,:), AP_NAME(I_MOMY), AP_DESC(I_MOMY), AP_UNIT(I_MOMY), ndim=3 )
    call MONIT_in( RHOT(:,:,:), AP_NAME(I_RHOT), AP_DESC(I_RHOT), AP_UNIT(I_RHOT), ndim=3 )
    do iq = 1, QA
       RHOQ(:,:,:) = DENS(:,:,:) * QTRC(:,:,:,iq)

       call MONIT_in( RHOQ(:,:,:), AQ_NAME(iq), AQ_DESC(iq), AQ_UNIT(iq), ndim=3 )
    enddo

    call THERMODYN_qd( QDRY(:,:,:),  & ! [OUT]
                       QTRC(:,:,:,:) ) ! [IN]

    call THERMODYN_temp_pres( TEMP(:,:,:),  & ! [OUT]
                              PRES(:,:,:),  & ! [OUT]
                              DENS(:,:,:),  & ! [IN]
                              RHOT(:,:,:),  & ! [IN]
                              QTRC(:,:,:,:) ) ! [IN]

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       VELZ(k,i,j) = 0.5_RP * ( MOMZ(k-1,i,j)+MOMZ(k,i,j) ) / DENS(k,i,j)
       VELX(k,i,j) = 0.5_RP * ( MOMX(k,i-1,j)+MOMX(k,i,j) ) / DENS(k,i,j)
       VELY(k,i,j) = 0.5_RP * ( MOMY(k,i,j-1)+MOMY(k,i,j) ) / DENS(k,i,j)

       ENGP(k,i,j) = DENS(k,i,j) * GRAV * CZ(k)

       ENGK(k,i,j) = 0.5_RP * DENS(k,i,j) * ( VELZ(k,i,j)**2 &
                                            + VELX(k,i,j)**2 &
                                            + VELY(k,i,j)**2 )

       ENGI(k,i,j) = DENS(k,i,j) * QDRY(k,i,j) * TEMP(k,i,j) * CVdry
       do iq = QQS, QQE
          ENGI(k,i,j) = ENGI(k,i,j) &
                      + DENS(k,i,j) * QTRC(k,i,j,iq) * TEMP(k,i,j) * CVw(iq)
       enddo

       ENGT(k,i,j) = ENGP(k,i,j) + ENGK(k,i,j) + ENGI(k,i,j)
    enddo
    enddo
    enddo

    RHOQ(:,:,:) = DENS(:,:,:) * QDRY (:,:,:)

    call MONIT_put( AD_MONIT_id(I_QDRY), RHOQ(:,:,:) )

    RHOQ(:,:,:) = DENS(:,:,:) * ( 1.0_RP - QDRY (:,:,:) )

    call MONIT_put( AD_MONIT_id(I_QTOT), RHOQ(:,:,:) )

    call MONIT_put( AD_MONIT_id(I_ENGT), ENGT(:,:,:) )
    call MONIT_put( AD_MONIT_id(I_ENGP), ENGP(:,:,:) )
    call MONIT_put( AD_MONIT_id(I_ENGK), ENGK(:,:,:) )
    call MONIT_put( AD_MONIT_id(I_ENGI), ENGI(:,:,:) )

    call MONIT_write('MAIN')

    return
  end subroutine ATMOS_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Write restart of atmospheric variables
  subroutine ATMOS_vars_restart_write
    use scale_process, only: &
       PRC_master, &
       PRC_myrank, &
       PRC_2Drank
    use scale_time, only: &
       NOWSEC => TIME_NOWDAYSEC
    use gtool_file_h, only: &
       File_REAL4, &
       File_REAL8
    use gtool_file, only: &
       FileCreate, &
       FileAddVariable, &
       FilePutAxis, &
       FilePutAssociatedCoordinates, &
       FileWrite, &
       FileClose
    use scale_grid, only: &
       GRID_CZ,    &
       GRID_CX,    &
       GRID_CY,    &
       GRID_FZ,    &
       GRID_FX,    &
       GRID_FY,    &
       GRID_CDZ,   &
       GRID_CDX,   &
       GRID_CDY,   &
       GRID_FDZ,   &
       GRID_FDX,   &
       GRID_FDY,   &
       GRID_CBFZ,  &
       GRID_CBFX,  &
       GRID_CBFY,  &
       GRID_FBFZ,  &
       GRID_FBFX,  &
       GRID_FBFY,  &
       GRID_CXG,   &
       GRID_CYG,   &
       GRID_FXG,   &
       GRID_FYG,   &
       GRID_CBFXG, &
       GRID_CBFYG, &
       GRID_FBFXG, &
       GRID_FBFYG
    use scale_grid_real, only: &
       REAL_LON, &
       REAL_LONX, &
       REAL_LAT, &
       REAL_LATY
#ifdef _SDM
    use scale_atmos_phy_mp_sdm, only: &
       sd_rest_flg_out, &
       ATMOS_PHY_MP_sdm_restart_out
#endif
    implicit none

    real(RP) :: restart_atmos(KMAX,IMAX,JMAX) !> restart file (no HALO)

    character(len=H_LONG) :: basename

    integer :: fid, ap_vid(5), aq_vid(QA)
    logical :: fileexisted
    integer :: dtype
    integer :: iq
    integer :: n

    integer :: rankidx(2)
    !---------------------------------------------------------------------------

#ifdef _SDM
    if( sd_rest_flg_out ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Output random number for SDM ***'
       call ATMOS_PHY_MP_sdm_restart_out(NOWSEC)
    endif
#endif

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Output restart file (atmos) ***'

    call PROF_rapstart('FILE O NetCDF')

    basename = ''
    write(basename(1:15), '(F15.3)') NOWSEC
    do n = 1, 15
       if ( basename(n:n) == ' ' ) basename(n:n) = '0'
    end do
    basename = trim(ATMOS_RESTART_OUT_BASENAME) // '_' // trim(basename)

    rankidx(1) = PRC_2Drank(PRC_myrank,1)
    rankidx(2) = PRC_2Drank(PRC_myrank,2)

    if ( RP == 8 ) then
       dtype = File_REAL8
    else if ( RP == 4 ) then
       dtype = File_REAL4
    endif

    call FileCreate( fid,                     & ! (out)
                     fileexisted,             & ! (out)
                     basename,                & ! (in)
                     ATMOS_RESTART_OUT_TITLE, & ! (in)
                     H_SOURCE,                & ! (in)
                     H_INSTITUTE,             & ! (in)
                     PRC_master,              & ! (in)
                     PRC_myrank,              & ! (in)
                     rankidx                  ) ! (in)

    call FilePutAxis( fid, 'z', 'Z', 'm', 'z', dtype, GRID_CZ(KS:KE) )
    call FilePutAxis( fid, 'x', 'X', 'm', 'x', dtype, GRID_CX(IS:IE) )
    call FilePutAxis( fid, 'y', 'Y', 'm', 'y', dtype, GRID_CY(JS:JE) )

    call FilePutAxis( fid, 'zh', 'Z (half level)', 'm', 'zh', dtype, GRID_FZ(KS:KE) )
    call FilePutAxis( fid, 'xh', 'X (half level)', 'm', 'xh', dtype, GRID_FX(IS:IE) )
    call FilePutAxis( fid, 'yh', 'Y (half level)', 'm', 'yh', dtype, GRID_FY(JS:JE) )

    call FilePutAxis( fid, 'CZ', 'Grid Center Position Z', 'm', 'CZ', dtype, GRID_CZ )
    call FilePutAxis( fid, 'CX', 'Grid Center Position X', 'm', 'CX', dtype, GRID_CX )
    call FilePutAxis( fid, 'CY', 'Grid Center Position Y', 'm', 'CY', dtype, GRID_CY )
    call FilePutAxis( fid, 'FZ', 'Grid Face Position Z',   'm', 'FZ', dtype, GRID_FZ )
    call FilePutAxis( fid, 'FX', 'Grid Face Position X',   'm', 'FX', dtype, GRID_FX )
    call FilePutAxis( fid, 'FY', 'Grid Face Position Y',   'm', 'FY', dtype, GRID_FY )

    call FilePutAxis( fid, 'CDZ', 'Grid Cell length Z', 'm', 'CZ',  dtype, GRID_CDZ )
    call FilePutAxis( fid, 'CDX', 'Grid Cell length X', 'm', 'CX',  dtype, GRID_CDX )
    call FilePutAxis( fid, 'CDY', 'Grid Cell length Y', 'm', 'CY',  dtype, GRID_CDY )
    call FilePutAxis( fid, 'FDZ', 'Grid distance Z',    'm', 'FDZ', dtype, GRID_FDZ )
    call FilePutAxis( fid, 'FDX', 'Grid distance X',    'm', 'FDX', dtype, GRID_FDX )
    call FilePutAxis( fid, 'FDY', 'Grid distance Y',    'm', 'FDY', dtype, GRID_FDY )

    call FilePutAxis( fid, 'CBFZ', 'Boundary factor Center Z', '1', 'CZ', dtype, GRID_CBFZ )
    call FilePutAxis( fid, 'CBFX', 'Boundary factor Center X', '1', 'CX', dtype, GRID_CBFX )
    call FilePutAxis( fid, 'CBFY', 'Boundary factor Center Y', '1', 'CY', dtype, GRID_CBFY )
    call FilePutAxis( fid, 'FBFZ', 'Boundary factor Face Z',   '1', 'CZ', dtype, GRID_FBFZ )
    call FilePutAxis( fid, 'FBFX', 'Boundary factor Face X',   '1', 'CX', dtype, GRID_FBFX )
    call FilePutAxis( fid, 'FBFY', 'Boundary factor Face Y',   '1', 'CY', dtype, GRID_FBFY )

    call FilePutAxis( fid, 'CXG', 'Grid Center Position X (global)', 'm', 'CXG', dtype, GRID_CXG )
    call FilePutAxis( fid, 'CYG', 'Grid Center Position Y (global)', 'm', 'CYG', dtype, GRID_CYG )
    call FilePutAxis( fid, 'FXG', 'Grid Face Position X (global)',   'm', 'FXG', dtype, GRID_FXG )
    call FilePutAxis( fid, 'FYG', 'Grid Face Position Y (global)',   'm', 'FYG', dtype, GRID_FYG )

    call FilePutAxis( fid, 'CBFXG', 'Boundary factor Center X (global)', '1', 'CXG', dtype, GRID_CBFXG )
    call FilePutAxis( fid, 'CBFYG', 'Boundary factor Center Y (global)', '1', 'CYG', dtype, GRID_CBFYG )
    call FilePutAxis( fid, 'FBFXG', 'Boundary factor Face X (global)',   '1', 'CXG', dtype, GRID_FBFXG )
    call FilePutAxis( fid, 'FBFYG', 'Boundary factor Face Y (global)',   '1', 'CYG', dtype, GRID_FBFYG )

    call FilePutAssociatedCoordinates( fid, &
         'lon', 'longitude', 'degrees_east', (/'x', 'y'/), dtype, REAL_LON(IS:IE,JS:JE) )
    call FilePutAssociatedCoordinates( fid, &
         'lonh', 'longitude (half level)', 'degrees_east', (/'xh', 'y '/), &
         dtype, REAL_LONX(IS:IE,JS:JE) )
    call FilePutAssociatedCoordinates( fid, &
         'lat', 'latitude', 'degrees_north', (/'x', 'y'/), dtype, REAL_LAT(IS:IE,JS:JE) )
    call FilePutAssociatedCoordinates( fid, &
         'lath', 'latitude (half level)', 'degrees_north', (/'x ', 'yh'/), &
         dtype, REAL_LATY(IS:IE,JS:JE) )

    call FileAddVariable( ap_vid(I_DENS),                                         & ! (out)
                          fid, AP_NAME(I_DENS), AP_DESC(I_DENS), AP_UNIT(I_DENS), & ! (in)
                          (/'z  ','lon','lat'/), dtype                            ) ! (in)
    call FileAddVariable( ap_vid(I_MOMZ),                                         & ! (out)
                          fid, AP_NAME(I_MOMZ), AP_DESC(I_MOMZ), AP_UNIT(I_MOMZ), & ! (in)
                          (/'zh ','lon','lat'/), dtype                            ) ! (in)
    call FileAddVariable( ap_vid(I_MOMX),                                         & ! (out)
                          fid, AP_NAME(I_MOMX), AP_DESC(I_MOMX), AP_UNIT(I_MOMX), & ! (in)
                          (/'z   ','lonh','lat '/), dtype                          ) ! (in)
    call FileAddVariable( ap_vid(I_MOMY),                                         & ! (out)
                          fid, AP_NAME(I_MOMY), AP_DESC(I_MOMY), AP_UNIT(I_MOMY), & ! (in)
                          (/'z   ','lon ','lath'/), dtype                          ) ! (in)
    call FileAddVariable( ap_vid(I_RHOT),                                         & ! (out)
                          fid, AP_NAME(I_RHOT), AP_DESC(I_RHOT), AP_UNIT(I_RHOT), & ! (in)
                          (/'z  ','lon','lat'/), dtype                            ) ! (in)
    do iq = 1, QA
       call FileAddVariable( aq_vid(iq),                                 & ! (out)
                             fid, AQ_NAME(iq), AQ_DESC(iq), AQ_UNIT(iq), & ! (in)
                             (/'z  ','lon','lat'/), dtype                             ) ! (in)
    end do

    restart_atmos(1:KMAX,1:IMAX,1:JMAX) = DENS(KS:KE,IS:IE,JS:JE)
    call FileWrite( ap_vid(I_DENS), restart_atmos(:,:,:), NOWSEC, NOWSEC )

    restart_atmos(1:KMAX,1:IMAX,1:JMAX) = MOMZ(KS:KE,IS:IE,JS:JE)
    call FileWrite( ap_vid(I_MOMZ), restart_atmos(:,:,:), NOWSEC, NOWSEC )

    restart_atmos(1:KMAX,1:IMAX,1:JMAX) = MOMX(KS:KE,IS:IE,JS:JE)
    call FileWrite( ap_vid(I_MOMX), restart_atmos(:,:,:), NOWSEC, NOWSEC )

    restart_atmos(1:KMAX,1:IMAX,1:JMAX) = MOMY(KS:KE,IS:IE,JS:JE)
    call FileWrite( ap_vid(I_MOMY), restart_atmos(:,:,:), NOWSEC, NOWSEC )

    restart_atmos(1:KMAX,1:IMAX,1:JMAX) = RHOT(KS:KE,IS:IE,JS:JE)
    call FileWrite( ap_vid(I_RHOT), restart_atmos(:,:,:), NOWSEC, NOWSEC )

    do iq = 1, QA
       restart_atmos(1:KMAX,1:IMAX,1:JMAX) = QTRC(KS:KE,IS:IE,JS:JE,iq)
       call FileWrite( aq_vid(iq), restart_atmos(:,:,:), NOWSEC, NOWSEC )
    enddo

    call FileClose( fid )

    call PROF_rapend  ('FILE O NetCDF')

    call ATMOS_vars_total

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
       LH0   => CONST_LH0,   &
       P00   => CONST_PRE00
    use scale_time, only: &
       TIME_DTSEC
    use scale_grid, only: &
       CZ   => GRID_CZ,   &
       CDZ  => GRID_CDZ,  &
       RCDZ => GRID_RCDZ, &
       RCDX => GRID_RCDX, &
       RCDY => GRID_RCDY
    use scale_history, only: &
       HIST_in
    use scale_monitor, only: &
       MONIT_put, &
       MONIT_in
    use scale_atmos_thermodyn, only: &
       THERMODYN_qd => ATMOS_THERMODYN_qd, &
       CPw => AQ_CP,                       &
       CVw => AQ_CV
    use scale_atmos_saturation, only: &
       SATURATION_dens2qsat_all => ATMOS_SATURATION_dens2qsat_all, &
       SATURATION_dens2qsat_liq => ATMOS_SATURATION_dens2qsat_liq, &
       SATURATION_dens2qsat_ice => ATMOS_SATURATION_dens2qsat_ice
    implicit none

    real(RP) :: VELZ  (KA,IA,JA) ! velocity w at cell center [m/s]
    real(RP) :: VELX  (KA,IA,JA) ! velocity u at cell center [m/s]
    real(RP) :: VELY  (KA,IA,JA) ! velocity v at cell center [m/s]
    real(RP) :: POTT  (KA,IA,JA) ! potential temperature [K]

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
    real(RP) :: PRES  (KA,IA,JA) ! pressure    [Pa]
    real(RP) :: TEMP  (KA,IA,JA) ! temperature [K]

    real(RP) :: POTL  (KA,IA,JA) ! liquid water potential temperature [K]
    real(RP) :: RH    (KA,IA,JA) ! relative humidity (liquid+ice)      [%]
    real(RP) :: RHL   (KA,IA,JA) ! relative humidity against to liquid [%]
    real(RP) :: RHI   (KA,IA,JA) ! relative humidity against to ice    [%]

    real(RP) :: VOR   (KA,IA,JA) ! vertical vorticity    [1/s]
    real(RP) :: DIV   (KA,IA,JA) ! divergence            [1/s]
    real(RP) :: HDIV  (KA,IA,JA) ! horizontal divergence [1/s]

    real(RP) :: ENGT  (KA,IA,JA) ! total     energy [J/m3]
    real(RP) :: ENGP  (KA,IA,JA) ! potential energy [J/m3]
    real(RP) :: ENGK  (KA,IA,JA) ! kinetic   energy [J/m3]
    real(RP) :: ENGI  (KA,IA,JA) ! internal  energy [J/m3]

    real(RP) :: RHOQ  (KA,IA,JA)
    real(RP) :: QSAT  (KA,IA,JA)
    real(RP) :: VELXH (KA,IA,JA)
    real(RP) :: VELYH (KA,IA,JA)

    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    ! value check for prognostic variables
    if ( ATMOS_VARS_CHECKRANGE ) then
       call VALCHECK( DENS(:,:,:),    0.0_RP,    2.0_RP, AP_NAME(I_DENS), __FILE__, __LINE__ )
       call VALCHECK( MOMZ(:,:,:), -200.0_RP,  200.0_RP, AP_NAME(I_MOMZ), __FILE__, __LINE__ )
       call VALCHECK( MOMX(:,:,:), -200.0_RP,  200.0_RP, AP_NAME(I_MOMX), __FILE__, __LINE__ )
       call VALCHECK( MOMY(:,:,:), -200.0_RP,  200.0_RP, AP_NAME(I_MOMY), __FILE__, __LINE__ )
       call VALCHECK( RHOT(:,:,:),    0.0_RP, 1000.0_RP, AP_NAME(I_RHOT), __FILE__, __LINE__ )
    endif

    ! history output of prognostic variables
    call HIST_in( DENS(:,:,:), AP_NAME(I_DENS), AP_DESC(I_DENS), AP_UNIT(I_DENS), TIME_DTSEC )
    call HIST_in( MOMZ(:,:,:), AP_NAME(I_MOMZ), AP_DESC(I_MOMZ), AP_UNIT(I_MOMZ), TIME_DTSEC )
    call HIST_in( MOMX(:,:,:), AP_NAME(I_MOMX), AP_DESC(I_MOMX), AP_UNIT(I_MOMX), TIME_DTSEC )
    call HIST_in( MOMY(:,:,:), AP_NAME(I_MOMY), AP_DESC(I_MOMY), AP_UNIT(I_MOMY), TIME_DTSEC )
    call HIST_in( RHOT(:,:,:), AP_NAME(I_RHOT), AP_DESC(I_RHOT), AP_UNIT(I_RHOT), TIME_DTSEC )
    do iq = 1, QA
       call HIST_in( QTRC(:,:,:,iq), AQ_NAME(iq), AQ_DESC(iq), AQ_UNIT(iq), TIME_DTSEC )
    enddo

    ! monitor output of prognostic variables
    call MONIT_in( DENS(:,:,:), AP_NAME(I_DENS), AP_DESC(I_DENS), AP_UNIT(I_DENS), ndim=3 )
    call MONIT_in( MOMZ(:,:,:), AP_NAME(I_MOMZ), AP_DESC(I_MOMZ), AP_UNIT(I_MOMZ), ndim=3 )
    call MONIT_in( MOMX(:,:,:), AP_NAME(I_MOMX), AP_DESC(I_MOMX), AP_UNIT(I_MOMX), ndim=3 )
    call MONIT_in( MOMY(:,:,:), AP_NAME(I_MOMY), AP_DESC(I_MOMY), AP_UNIT(I_MOMY), ndim=3 )
    call MONIT_in( RHOT(:,:,:), AP_NAME(I_RHOT), AP_DESC(I_RHOT), AP_UNIT(I_RHOT), ndim=3 )
    do iq = 1, QA
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          RHOQ(k,i,j) = DENS(k,i,j) * QTRC(k,i,j,iq)
       enddo
       enddo
       enddo

       call MONIT_in( RHOQ(:,:,:), AQ_NAME(iq), AQ_DESC(iq), AQ_UNIT(iq), ndim=3 )
    enddo

    ! prepare and history output of diagnostic variables

    if ( AD_PREP_sw(I_VELZ) > 0 ) then
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          VELZ(k,i,j) = 0.5_RP * ( MOMZ(k-1,i,j)+MOMZ(k,i,j) ) / DENS(k,i,j)
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_VELX) > 0 ) then
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          VELX(k,i,j) = 0.5_RP * ( MOMX(k,i-1,j)+MOMX(k,i,j) ) / DENS(k,i,j)
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_VELY) > 0 ) then
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          VELY(k,i,j) = 0.5_RP * ( MOMY(k,i,j-1)+MOMY(k,i,j) ) / DENS(k,i,j)
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_POTT) > 0 ) then
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          POTT(k,i,j) = RHOT(k,i,j) / DENS(k,i,j)
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_QDRY) > 0 ) then
       call THERMODYN_qd( QDRY(:,:,:),  & ! [OUT]
                          QTRC(:,:,:,:) ) ! [IN]
    endif

    if ( AD_PREP_sw(I_QTOT) > 0 ) then
       do j  = JS, JE
       do i  = IS, IE
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
       do j  = JS, JE
       do i  = IS, IE
          LWP(i,j) = 0.0_RP
          do k  = KS, KE
             LWP(i,j) = LWP(i,j) + QLIQ(k,i,j) * DENS(k,i,j) * CDZ(k) * 1.E3_RP ! [kg/m2->g/m2]
          enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_IWP) > 0 ) then
       do j  = JS, JE
       do i  = IS, IE
          IWP(i,j) = 0.0_RP
          do k  = KS, KE
             IWP(i,j) = IWP(i,j) + QICE(k,i,j) * DENS(k,i,j) * CDZ(k) * 1.E3_RP ! [kg/m2->g/m2]
          enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_RTOT) > 0 ) then
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j  = JS, JE
       do i  = IS, IE
       do k  = KS, KE
          RTOT (k,i,j) = Rdry * QDRY(k,i,j) + Rvap * QTRC(k,i,j,I_QV)
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_CPTOT) > 0 ) then
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j  = JS, JE
       do i  = IS, IE
       do k  = KS, KE
          CPTOT(k,i,j) = CPdry * QDRY(k,i,j)
       enddo
       enddo
       enddo

       do iq = QQS, QQE
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j  = JS, JE
       do i  = IS, IE
       do k  = KS, KE
          CPTOT(k,i,j) = CPTOT(k,i,j) + QTRC(k,i,j,iq) * CPw(iq)
       enddo
       enddo
       enddo
       enddo

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j  = JS, JE
       do i  = IS, IE
       do k  = KS, KE
          CPovCV(k,i,j) = CPTOT(k,i,j) / ( CPTOT(k,i,j) - RTOT(k,i,j) )
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_PRES) > 0 ) then
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          PRES(k,i,j) = P00 * ( RHOT(k,i,j) * RTOT(k,i,j) / P00 )**CPovCV(k,i,j)
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_TEMP) > 0 ) then
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          TEMP(k,i,j) = PRES(k,i,j) / ( DENS(k,i,j) * RTOT(k,i,j) )
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_POTL) > 0 ) then
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          POTL(k,i,j) = POTT(k,i,j) &
                      - LH0 / CPdry * QLIQ(k,i,j) * POTT(k,i,j) / TEMP(k,i,j)
       enddo
       enddo
       enddo
    endif

    if ( AD_HIST_id(I_RH) > 0 ) then
       call SATURATION_dens2qsat_all( QSAT(:,:,:), & ! [OUT]
                                      TEMP(:,:,:), & ! [IN]
                                      DENS(:,:,:)  ) ! [IN]

       do j  = JS, JE
       do i  = IS, IE
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
       do j  = JS, JE
       do i  = IS, IE
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
       do j  = JS, JE
       do i  = IS, IE
       do k  = KS, KE
          RHI(k,i,j) = QTRC(k,i,j,I_QV) / QSAT(k,i,j) * 1.0E2_RP
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_VOR) > 0 ) then
       ! at u, v, layer
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS-1, JE
       do i = IS-1, IE
       do k = KS, KE
          VELXH(k,i,j) = 2.0_RP * ( MOMX(k,i,j)+MOMX(k,i,j+1) )                             &
                                / ( DENS(k,i,j)+DENS(k,i+1,j)+DENS(k,i,j+1)+DENS(k,i+1,j+1) )
       enddo
       enddo
       enddo

       ! at u, v, layer
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS-1, JE
       do i = IS-1, IE
       do k = KS, KE
          VELYH(k,i,j) = 2.0_RP * ( MOMY(k,i,j)+MOMY(k,i+1,j) )                             &
                                / ( DENS(k,i,j)+DENS(k,i+1,j)+DENS(k,i,j+1)+DENS(k,i+1,j+1) )
       enddo
       enddo
       enddo

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          VOR(k,i,j) = 0.5_RP * ( ( VELYH(k,i,j  ) - VELYH(k,i-1,j  ) ) * RCDX(i) &
                                + ( VELYH(k,i,j-1) - VELYH(k,i-1,j-1) ) * RCDX(i) &
                                - ( VELXH(k,i  ,j) - VELXH(k,i  ,j-1) ) * RCDY(j) &
                                - ( VELXH(k,i-1,j) - VELXH(k,i-1,j-1) ) * RCDY(j) )
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_HDIV) > 0 ) then
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          HDIV(k,i,j) = ( MOMX(k,i,j) - MOMX(k  ,i-1,j  ) ) * RCDX(i) &
                      + ( MOMY(k,i,j) - MOMY(k  ,i  ,j-1) ) * RCDY(j)
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_DIV) > 0 ) then
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          DIV(k,i,j) = ( MOMZ(k,i,j) - MOMZ(k-1,i  ,j  ) ) * RCDZ(k) &
                     + ( MOMX(k,i,j) - MOMX(k  ,i-1,j  ) ) * RCDX(i) &
                     + ( MOMY(k,i,j) - MOMY(k  ,i  ,j-1) ) * RCDY(j)
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_ENGP) > 0 ) then
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          ENGP(k,i,j) = DENS(k,i,j) * GRAV * CZ(k)
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_ENGK) > 0 ) then
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          ENGK(k,i,j) = 0.5_RP * DENS(k,i,j) * ( VELZ(k,i,j)**2 &
                                               + VELX(k,i,j)**2 &
                                               + VELY(k,i,j)**2 )
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_ENGI) > 0 ) then
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          ENGI(k,i,j) = DENS(k,i,j) * QDRY(k,i,j) * TEMP(k,i,j) * CVdry
       enddo
       enddo
       enddo

       do iq = QQS, QQE
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(3)
       do j  = JS, JE
       do i  = IS, IE
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
       do j  = JS, JE
       do i  = IS, IE
       do k  = KS, KE
          ENGT(k,i,j) = ENGP(k,i,j) + ENGK(k,i,j) + ENGI(k,i,j)
       enddo
       enddo
       enddo
    endif

    call HIST_in( VELZ (:,:,:), 'W',     'velocity w',             'm/s',    TIME_DTSEC )
    call HIST_in( VELX (:,:,:), 'U',     'velocity u',             'm/s',    TIME_DTSEC )
    call HIST_in( VELY (:,:,:), 'V',     'velocity v',             'm/s',    TIME_DTSEC )
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

    call HIST_in( ENGT (:,:,:), 'ENGT',  'total energy',           'J/m3',   TIME_DTSEC )
    call HIST_in( ENGP (:,:,:), 'ENGP',  'potential energy',       'J/m3',   TIME_DTSEC )
    call HIST_in( ENGK (:,:,:), 'ENGK',  'kinetic energy',         'J/m3',   TIME_DTSEC )
    call HIST_in( ENGI (:,:,:), 'ENGI',  'internal energy',        'J/m3',   TIME_DTSEC )

    if ( AD_MONIT_id(I_QDRY) > 0 ) then
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          RHOQ(k,i,j) = DENS(k,i,j) * QDRY(k,i,j)
       enddo
       enddo
       enddo
       call MONIT_put( AD_MONIT_id(I_QDRY), RHOQ(:,:,:) )
    endif

    if ( AD_MONIT_id(I_QTOT) > 0 ) then
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          RHOQ(k,i,j) = DENS(k,i,j) * QTOT(k,i,j)
       enddo
       enddo
       enddo
       call MONIT_put( AD_MONIT_id(I_QTOT), RHOQ(:,:,:) )
    endif

    call MONIT_put( AD_MONIT_id(I_ENGT), ENGT(:,:,:) )
    call MONIT_put( AD_MONIT_id(I_ENGP), ENGP(:,:,:) )
    call MONIT_put( AD_MONIT_id(I_ENGK), ENGK(:,:,:) )
    call MONIT_put( AD_MONIT_id(I_ENGI), ENGI(:,:,:) )

    return
  end subroutine ATMOS_vars_history

  !-----------------------------------------------------------------------------
  !> Budget monitor for atmosphere
  subroutine ATMOS_vars_total
    use scale_const, only: &
       GRAV   => CONST_GRAV,   &
       CVdry  => CONST_CVdry
    use scale_grid, only: &
       CZ   => GRID_CZ
    use scale_stats, only: &
       STAT_checktotal, &
       STAT_total
    use scale_atmos_thermodyn, only: &
       THERMODYN_qd        => ATMOS_THERMODYN_qd,        &
       THERMODYN_temp_pres => ATMOS_THERMODYN_temp_pres, &
       CVw => AQ_CV
    implicit none

    real(RP) :: VELZ  (KA,IA,JA) ! velocity w at cell center [m/s]
    real(RP) :: VELX  (KA,IA,JA) ! velocity u at cell center [m/s]
    real(RP) :: VELY  (KA,IA,JA) ! velocity v at cell center [m/s]

    real(RP) :: QDRY  (KA,IA,JA) ! dry air     [kg/kg]
    real(RP) :: PRES  (KA,IA,JA) ! pressure    [Pa]
    real(RP) :: TEMP  (KA,IA,JA) ! temperature [K]

    real(RP) :: ENGT  (KA,IA,JA) ! total     energy [J/m3]
    real(RP) :: ENGP  (KA,IA,JA) ! potential energy [J/m3]
    real(RP) :: ENGK  (KA,IA,JA) ! kinetic   energy [J/m3]
    real(RP) :: ENGI  (KA,IA,JA) ! internal  energy [J/m3]

    real(RP) :: RHOQ  (KA,IA,JA)

    real(RP) :: total ! dummy
    integer :: i, j, k, iq
    !---------------------------------------------------------------------------

    if ( STAT_checktotal ) then

       call STAT_total( total, DENS(:,:,:), AP_NAME(I_DENS) )
       call STAT_total( total, MOMZ(:,:,:), AP_NAME(I_MOMZ) )
       call STAT_total( total, MOMX(:,:,:), AP_NAME(I_MOMX) )
       call STAT_total( total, MOMY(:,:,:), AP_NAME(I_MOMY) )
       call STAT_total( total, RHOT(:,:,:), AP_NAME(I_RHOT) )
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

       call STAT_total( total, RHOQ(:,:,:), 'Qdry    ' )

       RHOQ(:,:,:) = DENS(:,:,:) * ( 1.0_RP - QDRY (:,:,:) ) ! Qtotal

       call STAT_total( total, RHOQ(:,:,:), 'Qtotal  ' )

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          VELZ(k,i,j) = 0.5_RP * ( MOMZ(k-1,i,j)+MOMZ(k,i,j) ) / DENS(k,i,j)
          VELX(k,i,j) = 0.5_RP * ( MOMX(k,i-1,j)+MOMX(k,i,j) ) / DENS(k,i,j)
          VELY(k,i,j) = 0.5_RP * ( MOMY(k,i,j-1)+MOMY(k,i,j) ) / DENS(k,i,j)

          ENGP(k,i,j) = DENS(k,i,j) * GRAV * CZ(k)

          ENGK(k,i,j) = 0.5_RP * DENS(k,i,j) * ( VELZ(k,i,j)**2 &
                                               + VELX(k,i,j)**2 &
                                               + VELY(k,i,j)**2 )

          ENGI(k,i,j) = DENS(k,i,j) * QDRY(k,i,j) * TEMP(k,i,j) * CVdry
          do iq = QQS, QQE
             ENGI(k,i,j) = ENGI(k,i,j) &
                         + DENS(k,i,j) * QTRC(k,i,j,iq) * TEMP(k,i,j) * CVw(iq)
          enddo

          ENGT(k,i,j) = ENGP(k,i,j) + ENGK(k,i,j) + ENGI(k,i,j)
       enddo
       enddo
       enddo

       call STAT_total( total, ENGT(:,:,:), 'ENGT    ' )
       call STAT_total( total, ENGP(:,:,:), 'ENGP    ' )
       call STAT_total( total, ENGK(:,:,:), 'ENGK    ' )
       call STAT_total( total, ENGI(:,:,:), 'ENGI    ' )

    endif

    return
  end subroutine ATMOS_vars_total

end module mod_atmos_vars
