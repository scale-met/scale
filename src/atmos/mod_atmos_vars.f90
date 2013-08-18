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
  use mod_stdio, only: &
     IO_FID_LOG, &
     IO_L,       &
     IO_SYSCHR,  &
     IO_FILECHR
  use mod_time, only: &
     TIME_rapstart, &
     TIME_rapend
  use gtool_file_h, only: &
     File_HLONG
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
# include "scale-les.h"
  include 'inc_precision.h'
  include 'inc_index.h'
  include 'inc_tracer.h'

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  character(len=IO_SYSCHR), public, save :: ATMOS_TYPE_DYN    = 'NONE'
  character(len=IO_SYSCHR), public, save :: ATMOS_TYPE_PHY_SF = 'NONE'
  character(len=IO_SYSCHR), public, save :: ATMOS_TYPE_PHY_TB = 'NONE'
  character(len=IO_SYSCHR), public, save :: ATMOS_TYPE_PHY_MP = 'NONE'
  character(len=IO_SYSCHR), public, save :: ATMOS_TYPE_PHY_RD = 'NONE'

  logical,                  public, save :: ATMOS_sw_dyn
  logical,                  public, save :: ATMOS_sw_phy_sf
  logical,                  public, save :: ATMOS_sw_phy_tb
  logical,                  public, save :: ATMOS_sw_phy_mp
  logical,                  public, save :: ATMOS_sw_phy_rd
  logical,                  public, save :: ATMOS_sw_restart
  logical,                  public, save :: ATMOS_sw_check

  logical,                  public, save :: ATMOS_USE_AVERAGE = .false.

  ! prognostic variables
  real(RP), public, target  :: DENS(KA,IA,JA)    ! Density    [kg/m3]
  real(RP), public, target  :: MOMZ(KA,IA,JA)    ! momentum z [kg/s/m2]
  real(RP), public, target  :: MOMX(KA,IA,JA)    ! momentum x [kg/s/m2]
  real(RP), public, target  :: MOMY(KA,IA,JA)    ! momentum y [kg/s/m2]
  real(RP), public, target  :: RHOT(KA,IA,JA)    ! DENS * POTT [K*kg/m3]
  real(RP), public, target  :: QTRC(KA,IA,JA,QA) ! ratio of mass of tracer to total mass[kg/kg]

  real(RP), public, target  :: DENS_avw(KA,IA,JA)
  real(RP), public, target  :: MOMZ_avw(KA,IA,JA)
  real(RP), public, target  :: MOMX_avw(KA,IA,JA)
  real(RP), public, target  :: MOMY_avw(KA,IA,JA)
  real(RP), public, target  :: RHOT_avw(KA,IA,JA)
  real(RP), public, target  :: QTRC_avw(KA,IA,JA,QA)

  real(RP), public, pointer :: DENS_av(:,:,:)
  real(RP), public, pointer :: MOMZ_av(:,:,:)
  real(RP), public, pointer :: MOMX_av(:,:,:)
  real(RP), public, pointer :: MOMY_av(:,:,:)
  real(RP), public, pointer :: RHOT_av(:,:,:)
  real(RP), public, pointer :: QTRC_av(:,:,:,:)

  ! tendency by physical processes
  real(RP), public, save    :: DENS_tp(KA,IA,JA)
  real(RP), public, save    :: MOMZ_tp(KA,IA,JA)
  real(RP), public, save    :: MOMX_tp(KA,IA,JA)
  real(RP), public, save    :: MOMY_tp(KA,IA,JA)
  real(RP), public, save    :: RHOT_tp(KA,IA,JA)
  real(RP), public, save    :: QTRC_tp(KA,IA,JA,QA)

  integer, public, parameter :: ZDIR = 1
  integer, public, parameter :: XDIR = 2
  integer, public, parameter :: YDIR = 3

  integer, public, parameter :: I_DENS = 1
  integer, public, parameter :: I_MOMZ = 2
  integer, public, parameter :: I_MOMX = 3
  integer, public, parameter :: I_MOMY = 4
  integer, public, parameter :: I_RHOT = 5
  integer, public, parameter :: I_QTRC = 6

  character(len=16), public, save :: AP_NAME(5)
  character(len=64), public, save :: AP_DESC(5)
  character(len=16), public, save :: AP_UNIT(5)

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
  logical,                   private, save :: ATMOS_RESTART_OUTPUT           = .false.
  character(len=IO_FILECHR), private, save :: ATMOS_RESTART_IN_BASENAME      = 'restart_in'
  character(len=IO_FILECHR), private, save :: ATMOS_RESTART_OUT_BASENAME     = 'restart_out'
  character(len=File_HLONG), private, save :: ATMOS_RESTART_OUT_TITLE        = 'SCALE3 PROGNOSTIC VARS.'
  character(len=File_HLONG), private, save :: ATMOS_RESTART_OUT_SOURCE       = 'SCALE-LES ver. '//VERSION
  character(len=File_HLONG), private, save :: ATMOS_RESTART_OUT_INSTITUTE    = 'AICS/RIKEN'
  logical,                   private, save :: ATMOS_RESTART_IN_ALLOWMISSINGQ = .false.

  logical,                   private, save :: ATMOS_RESTART_CHECK            = .false.
  character(len=IO_FILECHR), private, save :: ATMOS_RESTART_CHECK_BASENAME   = 'restart_check'
  real(RP),                  private, save :: ATMOS_RESTART_CHECK_CRITERION  = 1.E-6_RP

  logical,                   private, save :: ATMOS_VARS_CHECKRANGE          = .false.

  ! history output of prognostic variables
  integer, private, save      :: AP_HIST_id(5)
  integer, private, save      :: AQ_HIST_id(QA)

  ! history & monitor output of diagnostic variables
  integer, private, parameter :: AD_nmax = 23 ! number of diagnostic variables for history output

  integer, private, parameter :: I_VELZ  =  1 ! velocity w at cell center
  integer, private, parameter :: I_VELX  =  2 ! velocity u at cell center
  integer, private, parameter :: I_VELY  =  3 ! velocity v at cell center
  integer, private, parameter :: I_POTT  =  4 ! potential temperature

  integer, private, parameter :: I_QDRY  =  5 ! ratio of dry air            to total mass
  integer, private, parameter :: I_QTOT  =  6 ! ratio of total tracer       to total mass
  integer, private, parameter :: I_QHYD  =  7 ! ratio of total hydrometeor  to total mass
  integer, private, parameter :: I_QLIQ  =  8 ! ratio of total liquid water to total mass
  integer, private, parameter :: I_QICE  =  9 ! ratio of total ice    water to total mass

  integer, private, parameter :: I_RTOT  = 10 ! total gas constant
  integer, private, parameter :: I_CPTOT = 11 ! total heat capacity (constant pressure)
  integer, private, parameter :: I_PRES  = 12 ! pressure
  integer, private, parameter :: I_TEMP  = 13 ! temperature

  integer, private, parameter :: I_RH    = 14 ! relative humidity (liquid+ice)
  integer, private, parameter :: I_RHL   = 15 ! relative humidity against to liquid
  integer, private, parameter :: I_RHI   = 16 ! relative humidity against to ice

  integer, private, parameter :: I_VOR   = 17 ! vertical vorticity
  integer, private, parameter :: I_DIV   = 18 ! divergence
  integer, private, parameter :: I_HDIV  = 19 ! horizontal divergence

  integer, private, parameter :: I_ENGP  = 20 ! potential energy
  integer, private, parameter :: I_ENGK  = 21 ! kinetic   energy
  integer, private, parameter :: I_ENGI  = 22 ! internal  energy
  integer, private, parameter :: I_ENGT  = 23 ! total     energy

  integer, private, save      :: AD_HIST_id (AD_nmax)
  integer, private, save      :: AD_PREP_sw (AD_nmax)
  integer, private, save      :: AD_MONIT_id(AD_nmax)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_vars_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    use mod_history, only: &
       HIST_reg
    use mod_monitor, only: &
       MONIT_reg
    implicit none

    NAMELIST / PARAM_ATMOS / &
       ATMOS_TYPE_DYN,    &
       ATMOS_TYPE_PHY_SF, &
       ATMOS_TYPE_PHY_TB, &
       ATMOS_TYPE_PHY_MP, &
       ATMOS_TYPE_PHY_RD, &
       ATMOS_TYPE_PHY_AE, &
       ATMOS_TYPE_PHY_CH, &
       ATMOS_USE_AVERAGE

    NAMELIST / PARAM_ATMOS_VARS / &
       ATMOS_RESTART_IN_BASENAME,      &
       ATMOS_RESTART_IN_ALLOWMISSINGQ, &
       ATMOS_RESTART_OUTPUT,           &
       ATMOS_RESTART_OUT_BASENAME,     &
       ATMOS_RESTART_OUT_TITLE,        &
       ATMOS_RESTART_OUT_SOURCE,       &
       ATMOS_RESTART_OUT_INSTITUTE,    &
       ATMOS_RESTART_CHECK,            &
       ATMOS_RESTART_CHECK_BASENAME,   &
       ATMOS_RESTART_CHECK_CRITERION,  &
       ATMOS_VARS_CHECKRANGE

    integer :: ierr
    integer :: ip, iq
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Variables]/Categ[ATMOS]'

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

    if ( ATMOS_TYPE_DYN /= 'OFF' .AND. ATMOS_TYPE_DYN /= 'NONE' ) then
       if( IO_L ) write(IO_FID_LOG,*) '  Dynamical core   : ON'
       if( IO_L ) write(IO_FID_LOG,*) '  Tracer advection : ON'
       ATMOS_sw_dyn = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '  Dynamical core   : OFF'
       if( IO_L ) write(IO_FID_LOG,*) '  Tracer advection : OFF'
       ATMOS_sw_dyn = .false.
    endif

    if( IO_L ) write(IO_FID_LOG,*) 'Physics...'

    if ( ATMOS_TYPE_PHY_SF /= 'OFF' .AND. ATMOS_TYPE_PHY_SF /= 'NONE' ) then
       if( IO_L ) write(IO_FID_LOG,*) '  Surface Flux : ON'
       ATMOS_sw_phy_sf = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '  Surface Flux : OFF'
       ATMOS_sw_phy_sf = .false.
    endif
    if ( ATMOS_TYPE_PHY_TB /= 'OFF' .AND. ATMOS_TYPE_PHY_TB /= 'NONE' ) then
       if( IO_L ) write(IO_FID_LOG,*) '  Sub-grid Turbulence : ON'
       ATMOS_sw_phy_tb = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '  Sub-grid Turbulence : OFF'
       ATMOS_sw_phy_tb = .false.
    endif
    if ( ATMOS_TYPE_PHY_MP /= 'OFF' .AND. ATMOS_TYPE_PHY_MP /= 'NONE' ) then
       if( IO_L ) write(IO_FID_LOG,*) '  Cloud Microphysics  : ON'
       ATMOS_sw_phy_mp = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '  Cloud Microphysics  : OFF'
       ATMOS_sw_phy_mp = .false.
    endif
    if ( ATMOS_TYPE_PHY_RD /= 'OFF' .AND. ATMOS_TYPE_PHY_RD /= 'NONE' ) then
       if( IO_L ) write(IO_FID_LOG,*) '  Radiative transfer  : ON'
       ATMOS_sw_phy_rd = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '  Radiative transfer  : OFF'
       ATMOS_sw_phy_rd = .false.
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

    call HIST_reg( AP_HIST_id(I_DENS), AP_NAME(I_DENS), AP_DESC(I_DENS), AP_UNIT(I_DENS), ndim=3 )
    call HIST_reg( AP_HIST_id(I_MOMZ), AP_NAME(I_MOMZ), AP_DESC(I_MOMZ), AP_UNIT(I_MOMZ), ndim=3, zdim='half' )
    call HIST_reg( AP_HIST_id(I_MOMX), AP_NAME(I_MOMX), AP_DESC(I_MOMX), AP_UNIT(I_MOMX), ndim=3, xdim='half' )
    call HIST_reg( AP_HIST_id(I_MOMY), AP_NAME(I_MOMY), AP_DESC(I_MOMY), AP_UNIT(I_MOMY), ndim=3, ydim='half' )
    call HIST_reg( AP_HIST_id(I_RHOT), AP_NAME(I_RHOT), AP_DESC(I_RHOT), AP_UNIT(I_RHOT), ndim=3 )
    do iq = 1, QA
       call HIST_reg( AQ_HIST_id(iq), AQ_NAME(iq), AQ_DESC(iq), AQ_UNIT(iq), ndim=3 )
    enddo

    call HIST_reg( AD_HIST_id(I_VELZ),  'W',     'velocity w',             'm/s',    ndim=3 )
    call HIST_reg( AD_HIST_id(I_VELX),  'U',     'velocity u',             'm/s',    ndim=3 )
    call HIST_reg( AD_HIST_id(I_VELY),  'V',     'velocity v',             'm/s',    ndim=3 )
    call HIST_reg( AD_HIST_id(I_POTT),  'PT',    'potential temp.',        'K',      ndim=3 )

    call HIST_reg( AD_HIST_id(I_QDRY),  'QDRY',  'dry air',                'kg/kg',  ndim=3 )
    call HIST_reg( AD_HIST_id(I_QTOT),  'QTOT',  'total water',            'kg/kg',  ndim=3 )
    call HIST_reg( AD_HIST_id(I_QHYD),  'QHYD',  'total hydrometeors',     'kg/kg',  ndim=3 )
    call HIST_reg( AD_HIST_id(I_QLIQ),  'QLIQ',  'total liquid water',     'kg/kg',  ndim=3 )
    call HIST_reg( AD_HIST_id(I_QICE),  'QICE',  'total ice water',        'kg/kg',  ndim=3 )

    call HIST_reg( AD_HIST_id(I_RTOT),  'RTOT',  'Total gas constant',     'J/kg/K', ndim=3 )
    call HIST_reg( AD_HIST_id(I_CPTOT), 'CPTOT', 'Total heat capacity',    'J/kg/K', ndim=3 )
    call HIST_reg( AD_HIST_id(I_PRES),  'PRES',  'pressure',               'Pa',     ndim=3 )
    call HIST_reg( AD_HIST_id(I_TEMP),  'T',     'temperature',            'K',      ndim=3 )

    call HIST_reg( AD_HIST_id(I_RH),    'RH',    'relative humidity',      '%',      ndim=3 )
    call HIST_reg( AD_HIST_id(I_RHL),   'RHL',   'relative humidity(liq)', '%',      ndim=3 )
    call HIST_reg( AD_HIST_id(I_RHI),   'RHI',   'relative humidity(ice)', '%',      ndim=3 )

    call HIST_reg( AD_HIST_id(I_VOR),   'VOR',   'vertical vorticity',     '1/s',    ndim=3 )
    call HIST_reg( AD_HIST_id(I_DIV),   'DIV',   'divergence',             '1/s',    ndim=3 )
    call HIST_reg( AD_HIST_id(I_HDIV),  'HDIV',  'horizontal divergence',  '1/s',    ndim=3 )

    call HIST_reg( AD_HIST_id(I_ENGP),  'ENGP',  'potential energy',       'J/m3',   ndim=3 )
    call HIST_reg( AD_HIST_id(I_ENGK),  'ENGK',  'kinetic energy',         'J/m3',   ndim=3 )
    call HIST_reg( AD_HIST_id(I_ENGI),  'ENGI',  'internal energy',        'J/m3',   ndim=3 )

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
    use mod_comm, only: &
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
    use mod_process, only: &
       PRC_myrank
    use mod_const, only: &
       GRAV  => CONST_GRAV,  &
       CVdry => CONST_CVdry
    use mod_grid, only: &
       CZ => GRID_CZ
    use mod_monitor, only: &
       MONIT_put,  &
       MONIT_in,   &
       MONIT_write
    use mod_atmos_thermodyn, only: &
       THERMODYN_qd        => ATMOS_THERMODYN_qd,        &
       THERMODYN_temp_pres => ATMOS_THERMODYN_temp_pres, &
       CVw => AQ_CV
    implicit none

    real(RP) :: restart_atmos(KMAX,IMAX,JMAX) !> restart file (no HALO)

    character(len=IO_FILECHR) :: bname

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

    call TIME_rapstart('FILE I NetCDF')

    bname = ATMOS_RESTART_IN_BASENAME

    call FileRead( restart_atmos(:,:,:), bname, 'DENS', 1, PRC_myrank )
    DENS(KS:KE,IS:IE,JS:JE) = restart_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FileRead( restart_atmos(:,:,:), bname, 'MOMZ', 1, PRC_myrank )
    MOMZ(KS:KE,IS:IE,JS:JE) = restart_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FileRead( restart_atmos(:,:,:), bname, 'MOMX', 1, PRC_myrank )
    MOMX(KS:KE,IS:IE,JS:JE) = restart_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FileRead( restart_atmos(:,:,:), bname, 'MOMY', 1, PRC_myrank )
    MOMY(KS:KE,IS:IE,JS:JE) = restart_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FileRead( restart_atmos(:,:,:), bname, 'RHOT', 1, PRC_myrank )
    RHOT(KS:KE,IS:IE,JS:JE) = restart_atmos(1:KMAX,1:IMAX,1:JMAX)

    do iq = 1, QA
       call FileRead( restart_atmos(:,:,:), bname, AQ_NAME(iq), 1, PRC_myrank )
       QTRC(KS:KE,IS:IE,JS:JE,iq) = restart_atmos(1:KMAX,1:IMAX,1:JMAX)
    enddo

    call TIME_rapend  ('FILE I NetCDF')

    call ATMOS_vars_fillhalo

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
    use mod_process, only: &
       PRC_master, &
       PRC_myrank, &
       PRC_2Drank
    use mod_time, only: &
       NOWSEC => TIME_NOWDAYSEC
    use gtool_file_h, only: &
       File_REAL4, &
       File_REAL8
    use gtool_file, only: &
       FileCreate, &
       FileAddVariable, &
       FilePutAxis, &
       FilePutAdditionalAxis, &
       FileWrite, &
       FileClose
    use mod_grid, only: &
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
    implicit none

    real(RP) :: restart_atmos(KMAX,IMAX,JMAX) !> restart file (no HALO)

    character(len=IO_FILECHR) :: bname

    integer :: fid, ap_vid(5), aq_vid(QA)
    logical :: fileexisted
    integer :: dtype
    integer :: iq
    integer :: n

    integer :: rankidx(2)
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Output restart file (atmos) ***'

    call TIME_rapstart('FILE O NetCDF')

    bname = ''
    write(bname(1:15), '(F15.3)') NOWSEC
    do n = 1, 15
       if ( bname(n:n) == ' ' ) bname(n:n) = '0'
    end do
    bname = trim(ATMOS_RESTART_OUT_BASENAME) // '_' // trim(bname)

    if ( RP == 8 ) then
       dtype = File_REAL8
    else if ( RP == 4 ) then
       dtype = File_REAL4
    endif

    rankidx(1) = PRC_2Drank(PRC_myrank,1)
    rankidx(2) = PRC_2Drank(PRC_myrank,2)
    call FileCreate( fid, fileexisted,                          & ! (out)
         bname,                                                 & ! (in)
         ATMOS_RESTART_OUT_TITLE,                               & ! (in)
         ATMOS_RESTART_OUT_SOURCE,                              & ! (in)
         ATMOS_RESTART_OUT_INSTITUTE,                           & ! (in)
         (/'z','x','y'/), (/KMAX,IMAX,JMAX/), (/'Z','X','Y'/),  & ! (in)
         (/'m','m','m'/), (/dtype,dtype,dtype/),                & ! (in)
         PRC_master, PRC_myrank, rankidx                        ) ! (in)

    call FilePutAxis(fid, 'z', GRID_CZ(KS:KE))
    call FilePutAxis(fid, 'x', GRID_CX(IS:IE))
    call FilePutAxis(fid, 'y', GRID_CY(JS:JE))

    call FilePutAdditionalAxis(fid, 'zh', 'Z (half level)', 'm', 'zh', dtype, GRID_FZ(KS:KE))
    call FilePutAdditionalAxis(fid, 'xh', 'X (half level)', 'm', 'xh', dtype, GRID_FX(IS:IE))
    call FilePutAdditionalAxis(fid, 'yh', 'Y (half level)', 'm', 'yh', dtype, GRID_FY(JS:JE))

    call FilePutAdditionalAxis(fid, 'CZ', 'Grid Center Position Z', 'm', 'CZ', dtype, GRID_CZ)
    call FilePutAdditionalAxis(fid, 'CX', 'Grid Center Position X', 'm', 'CX', dtype, GRID_CX)
    call FilePutAdditionalAxis(fid, 'CY', 'Grid Center Position Y', 'm', 'CY', dtype, GRID_CY)
    call FilePutAdditionalAxis(fid, 'FZ', 'Grid Face Position Z', 'm', 'FZ', dtype, GRID_FZ)
    call FilePutAdditionalAxis(fid, 'FX', 'Grid Face Position X', 'm', 'FX', dtype, GRID_FX)
    call FilePutAdditionalAxis(fid, 'FY', 'Grid Face Position Y', 'm', 'FY', dtype, GRID_FY)

    call FilePutAdditionalAxis(fid, 'CDZ', 'Grid Cell length Z', 'm', 'CZ', dtype, GRID_CDZ)
    call FilePutAdditionalAxis(fid, 'CDX', 'Grid Cell length X', 'm', 'CX', dtype, GRID_CDX)
    call FilePutAdditionalAxis(fid, 'CDY', 'Grid Cell length Y', 'm', 'CY', dtype, GRID_CDY)
    call FilePutAdditionalAxis(fid, 'FDZ', 'Grid distance Z', 'm', 'FDZ', dtype, GRID_FDZ)
    call FilePutAdditionalAxis(fid, 'FDX', 'Grid distance X', 'm', 'FDX', dtype, GRID_FDX)
    call FilePutAdditionalAxis(fid, 'FDY', 'Grid distance Y', 'm', 'FDY', dtype, GRID_FDY)

    call FilePutAdditionalAxis(fid, 'CBFZ', 'Boundary factor Center Z', '1', 'CZ', dtype, GRID_CBFZ)
    call FilePutAdditionalAxis(fid, 'CBFX', 'Boundary factor Center X', '1', 'CX', dtype, GRID_CBFX)
    call FilePutAdditionalAxis(fid, 'CBFY', 'Boundary factor Center Y', '1', 'CY', dtype, GRID_CBFY)
    call FilePutAdditionalAxis(fid, 'FBFZ', 'Boundary factor Face Z', '1', 'CZ', dtype, GRID_FBFZ)
    call FilePutAdditionalAxis(fid, 'FBFX', 'Boundary factor Face X', '1', 'CX', dtype, GRID_FBFX)
    call FilePutAdditionalAxis(fid, 'FBFY', 'Boundary factor Face Y', '1', 'CY', dtype, GRID_FBFY)

    call FilePutAdditionalAxis(fid, 'CXG', 'Grid Center Position X (global)', 'm', 'CXG', dtype, GRID_CXG)
    call FilePutAdditionalAxis(fid, 'CYG', 'Grid Center Position Y (global)', 'm', 'CYG', dtype, GRID_CYG)
    call FilePutAdditionalAxis(fid, 'FXG', 'Grid Face Position X (global)', 'm', 'FXG', dtype, GRID_FXG)
    call FilePutAdditionalAxis(fid, 'FYG', 'Grid Face Position Y (global)', 'm', 'FYG', dtype, GRID_FYG)

    call FilePutAdditionalAxis(fid, 'CBFXG', 'Boundary factor Center X (global)', '1', 'CXG', dtype, GRID_CBFXG)
    call FilePutAdditionalAxis(fid, 'CBFYG', 'Boundary factor Center Y (global)', '1', 'CYG', dtype, GRID_CBFYG)
    call FilePutAdditionalAxis(fid, 'FBFXG', 'Boundary factor Face X (global)', '1', 'CXG', dtype, GRID_FBFXG)
    call FilePutAdditionalAxis(fid, 'FBFYG', 'Boundary factor Face Y (global)', '1', 'CYG', dtype, GRID_FBFYG)

    call FileAddVariable( ap_vid(I_DENS),                        & ! (out)
         fid, AP_NAME(I_DENS), AP_DESC(I_DENS), AP_UNIT(I_DENS), & ! (in)
         (/'z','x','y'/), dtype                                  ) ! (in)
    call FileAddVariable( ap_vid(I_MOMZ),                        & ! (out)
         fid, AP_NAME(I_MOMZ), AP_DESC(I_MOMZ), AP_UNIT(I_MOMZ), & ! (in)
         (/'zh','x ','y '/), dtype                               ) ! (in)
    call FileAddVariable( ap_vid(I_MOMX),                        & ! (out)
         fid, AP_NAME(I_MOMX), AP_DESC(I_MOMX), AP_UNIT(I_MOMX), & ! (in)
         (/'z ','xh','y '/), dtype                               ) ! (in)
    call FileAddVariable( ap_vid(I_MOMY),                        & ! (out)
         fid, AP_NAME(I_MOMY), AP_DESC(I_MOMY), AP_UNIT(I_MOMY), & ! (in)
         (/'z ','x ','yh'/), dtype                               ) ! (in)
    call FileAddVariable( ap_vid(I_RHOT),                        & ! (out)
         fid, AP_NAME(I_RHOT), AP_DESC(I_RHOT), AP_UNIT(I_RHOT), & ! (in)
         (/'z','x','y'/), dtype                                  ) ! (in)
    do iq = 1, QA
       call FileAddVariable( aq_vid(iq),                & ! (out)
            fid, AQ_NAME(iq), AQ_DESC(iq), AQ_UNIT(iq), & ! (in)
            (/'z','x','y'/), dtype                      ) ! (in)
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

    call TIME_rapend  ('FILE O NetCDF')

    call ATMOS_vars_total

    return
  end subroutine ATMOS_vars_restart_write

  !-----------------------------------------------------------------------------
  !> Check and compare between last data and sample data
  subroutine ATMOS_vars_restart_check
    use mod_process, only: &
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

    character(len=IO_FILECHR) :: bname

    logical :: datacheck
    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    call TIME_rapstart('Debug')

    write(*,*) 'Compare last Data with ', trim(ATMOS_RESTART_CHECK_BASENAME), 'on PE=', PRC_myrank
    write(*,*) '*** criterion = ', ATMOS_RESTART_CHECK_CRITERION
    datacheck = .true.

    bname = ATMOS_RESTART_CHECK_BASENAME

    call FileRead( restart_atmos(:,:,:), bname, 'DENS', 1, PRC_myrank )
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

    call FileRead( restart_atmos(:,:,:), bname, 'MOMZ', 1, PRC_myrank )
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

    call FileRead( restart_atmos(:,:,:), bname, 'MOMX', 1, PRC_myrank )
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

    call FileRead( restart_atmos(:,:,:), bname, 'MOMY', 1, PRC_myrank )
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

    call FileRead( restart_atmos(:,:,:), bname, 'RHOT', 1, PRC_myrank )
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
       call FileRead( restart_atmos(:,:,:), bname, AQ_NAME(iq), 1, PRC_myrank )
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

    call TIME_rapend('Debug')

    return
  end subroutine ATMOS_vars_restart_check

  !-----------------------------------------------------------------------------
  !> History output set for atmospheric variables
  subroutine ATMOS_vars_history
    use mod_const, only: &
       GRAV  => CONST_GRAV,  &
       Rdry  => CONST_Rdry,  &
       Rvap  => CONST_Rvap,  &
       CPdry => CONST_CPdry, &
       CVdry => CONST_CVdry, &
       P00   => CONST_PRE00
    use mod_misc, only: &
       MISC_valcheck
    use mod_time, only: &
       TIME_DTSEC
    use mod_grid, only: &
       CZ   => GRID_CZ,   &
       RCDZ => GRID_RCDZ, &
       RCDX => GRID_RCDX, &
       RCDY => GRID_RCDY
    use mod_comm, only: &
       COMM_horizontal_mean
    use mod_history, only: &
       HIST_put
    use mod_monitor, only: &
       MONIT_put, &
       MONIT_in
    use mod_atmos_thermodyn, only: &
       THERMODYN_qd => ATMOS_THERMODYN_qd, &
       CPw => AQ_CP,                       &
       CVw => AQ_CV
    use mod_atmos_saturation, only: &
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

    real(RP) :: RTOT  (KA,IA,JA) ! Total gas constant  [J/kg/K]
    real(RP) :: CPTOT (KA,IA,JA) ! Total heat capacity [J/kg/K]
    real(RP) :: CPovCV(KA,IA,JA) ! Cp/Cv
    real(RP) :: PRES  (KA,IA,JA) ! pressure    [Pa]
    real(RP) :: TEMP  (KA,IA,JA) ! temperature [K]

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
       call MISC_valcheck( DENS(:,:,:),    0.0_RP,    2.0_RP, AP_NAME(I_DENS) )
       call MISC_valcheck( MOMZ(:,:,:), -200.0_RP,  200.0_RP, AP_NAME(I_MOMZ) )
       call MISC_valcheck( MOMX(:,:,:), -200.0_RP,  200.0_RP, AP_NAME(I_MOMX) )
       call MISC_valcheck( MOMY(:,:,:), -200.0_RP,  200.0_RP, AP_NAME(I_MOMY) )
       call MISC_valcheck( RHOT(:,:,:),    0.0_RP, 1000.0_RP, AP_NAME(I_RHOT) )
    endif

    ! history output of prognostic variables
    call HIST_put( AP_HIST_id(I_DENS), DENS(:,:,:), TIME_DTSEC )
    call HIST_put( AP_HIST_id(I_MOMZ), MOMZ(:,:,:), TIME_DTSEC )
    call HIST_put( AP_HIST_id(I_MOMX), MOMX(:,:,:), TIME_DTSEC )
    call HIST_put( AP_HIST_id(I_MOMY), MOMY(:,:,:), TIME_DTSEC )
    call HIST_put( AP_HIST_id(I_RHOT), RHOT(:,:,:), TIME_DTSEC )
    do iq = 1, QA
       call HIST_put( AQ_HIST_id(iq), QTRC(:,:,:,iq), TIME_DTSEC )
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
          QTOT(k,i,j) = 1.D0 - QDRY(k,i,j)
       enddo
       enddo
       enddo
    endif

    if ( AD_PREP_sw(I_QHYD) > 0 ) then
       QHYD(:,:,:) = 0.0_RP
       if( I_QC > 0 ) QHYD(:,:,:) = QHYD(:,:,:) + QTRC(:,:,:,I_QC)
       if( I_QR > 0 ) QHYD(:,:,:) = QHYD(:,:,:) + QTRC(:,:,:,I_QR)
       if( I_QI > 0 ) QHYD(:,:,:) = QHYD(:,:,:) + QTRC(:,:,:,I_QI)
       if( I_QS > 0 ) QHYD(:,:,:) = QHYD(:,:,:) + QTRC(:,:,:,I_QS)
       if( I_QG > 0 ) QHYD(:,:,:) = QHYD(:,:,:) + QTRC(:,:,:,I_QG)
    endif

    if ( AD_PREP_sw(I_QLIQ) > 0 ) then
       QLIQ(:,:,:) = 0.0_RP
       if( I_QC > 0 ) QLIQ(:,:,:) = QLIQ(:,:,:) + QTRC(:,:,:,I_QC)
       if( I_QR > 0 ) QLIQ(:,:,:) = QLIQ(:,:,:) + QTRC(:,:,:,I_QR)
    endif

    if ( AD_PREP_sw(I_QICE) > 0 ) then
       QICE(:,:,:) = 0.0_RP
       if( I_QI > 0 ) QICE(:,:,:) = QICE(:,:,:) + QTRC(:,:,:,I_QI)
       if( I_QS > 0 ) QICE(:,:,:) = QICE(:,:,:) + QTRC(:,:,:,I_QS)
       if( I_QG > 0 ) QICE(:,:,:) = QICE(:,:,:) + QTRC(:,:,:,I_QG)
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
          RHI(k,i,j) = QTRC(k,i,j,I_QV) / QSAT(k,i,j) * 1.0E2_RP
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

    call HIST_put( AD_HIST_id(I_VELZ),  VELZ(:,:,:),  TIME_DTSEC )
    call HIST_put( AD_HIST_id(I_VELX),  VELX(:,:,:),  TIME_DTSEC )
    call HIST_put( AD_HIST_id(I_VELY),  VELY(:,:,:),  TIME_DTSEC )
    call HIST_put( AD_HIST_id(I_POTT),  POTT(:,:,:),  TIME_DTSEC )

    call HIST_put( AD_HIST_id(I_QDRY),  QDRY(:,:,:),  TIME_DTSEC )
    call HIST_put( AD_HIST_id(I_QTOT),  QTOT(:,:,:),  TIME_DTSEC )
    call HIST_put( AD_HIST_id(I_QHYD),  QHYD(:,:,:),  TIME_DTSEC )
    call HIST_put( AD_HIST_id(I_QLIQ),  QLIQ(:,:,:),  TIME_DTSEC )
    call HIST_put( AD_HIST_id(I_QICE),  QICE(:,:,:),  TIME_DTSEC )

    call HIST_put( AD_HIST_id(I_RTOT),  RTOT(:,:,:),  TIME_DTSEC )
    call HIST_put( AD_HIST_id(I_CPTOT), CPTOT(:,:,:), TIME_DTSEC )
    call HIST_put( AD_HIST_id(I_PRES),  PRES(:,:,:),  TIME_DTSEC )
    call HIST_put( AD_HIST_id(I_TEMP),  TEMP(:,:,:),  TIME_DTSEC )

    call HIST_put( AD_HIST_id(I_RH),    RH  (:,:,:),  TIME_DTSEC )
    call HIST_put( AD_HIST_id(I_RHL),   RHL (:,:,:),  TIME_DTSEC )
    call HIST_put( AD_HIST_id(I_RHI),   RHI (:,:,:),  TIME_DTSEC )

    call HIST_put( AD_HIST_id(I_VOR),   VOR (:,:,:),  TIME_DTSEC )
    call HIST_put( AD_HIST_id(I_DIV),   DIV (:,:,:),  TIME_DTSEC )
    call HIST_put( AD_HIST_id(I_HDIV),  HDIV(:,:,:),  TIME_DTSEC )

    call HIST_put( AD_HIST_id(I_ENGT),  ENGT(:,:,:),  TIME_DTSEC )
    call HIST_put( AD_HIST_id(I_ENGP),  ENGP(:,:,:),  TIME_DTSEC )
    call HIST_put( AD_HIST_id(I_ENGK),  ENGK(:,:,:),  TIME_DTSEC )
    call HIST_put( AD_HIST_id(I_ENGI),  ENGI(:,:,:),  TIME_DTSEC )

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
    use mod_const, only: &
       GRAV   => CONST_GRAV,   &
       CVdry  => CONST_CVdry
    use mod_grid, only: &
       CZ   => GRID_CZ
    use mod_comm, only: &
       COMM_total_doreport, &
       COMM_total
    use mod_atmos_thermodyn, only: &
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

    if ( COMM_total_doreport ) then

       call COMM_total( total, DENS(:,:,:), AP_NAME(I_DENS) )
       call COMM_total( total, MOMZ(:,:,:), AP_NAME(I_MOMZ) )
       call COMM_total( total, MOMX(:,:,:), AP_NAME(I_MOMX) )
       call COMM_total( total, MOMY(:,:,:), AP_NAME(I_MOMY) )
       call COMM_total( total, RHOT(:,:,:), AP_NAME(I_RHOT) )
       do iq = 1, QA
          RHOQ(:,:,:) = DENS(:,:,:) * QTRC(:,:,:,iq)

          call COMM_total( total, RHOQ(:,:,:), AQ_NAME(iq) )
       enddo

       call THERMODYN_qd( QDRY(:,:,:),  & ! [OUT]
                          QTRC(:,:,:,:) ) ! [IN]

       call THERMODYN_temp_pres( TEMP(:,:,:),  & ! [OUT]
                                 PRES(:,:,:),  & ! [OUT]
                                 DENS(:,:,:),  & ! [IN]
                                 RHOT(:,:,:),  & ! [IN]
                                 QTRC(:,:,:,:) ) ! [IN]

       RHOQ(:,:,:) = DENS(:,:,:) * QDRY (:,:,:)

       call COMM_total( total, RHOQ(:,:,:), 'Qdry    ' )

       RHOQ(:,:,:) = DENS(:,:,:) * ( 1.0_RP - QDRY (:,:,:) ) ! Qtotal

       call COMM_total( total, RHOQ(:,:,:), 'Qtotal  ' )

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

       call COMM_total( total, ENGT(:,:,:), 'ENGT    ' )
       call COMM_total( total, ENGP(:,:,:), 'ENGP    ' )
       call COMM_total( total, ENGK(:,:,:), 'ENGK    ' )
       call COMM_total( total, ENGI(:,:,:), 'ENGI    ' )

    endif

    return
  end subroutine ATMOS_vars_total

end module mod_atmos_vars
