!-------------------------------------------------------------------------------
!> module Atmospheric Variables
!!
!! @par Description
!!          Container for atmospheric variables
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2011-11-11 (H.Yashiro)  [new]
!! @li      2012-03-23 (H.Yashiro)  [mod] Explicit index parameter inclusion
!!
!<
!-------------------------------------------------------------------------------
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
  use mod_fileio_h, only: &
     FIO_HSHORT, &
     FIO_HMID
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
  include 'inc_index.h'
  include 'inc_tracer.h'

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(8), public, save :: DENS(KA,IA,JA)    ! Density    [kg/m3]
  real(8), public, save :: MOMZ(KA,IA,JA)    ! momentum z [kg/s/m2]
  real(8), public, save :: MOMX(KA,IA,JA)    ! momentum x [kg/s/m2]
  real(8), public, save :: MOMY(KA,IA,JA)    ! momentum y [kg/s/m2]
  real(8), public, save :: RHOT(KA,IA,JA)    ! DENS * POTT [K*kg/m3]
  real(8), public, save :: QTRC(KA,IA,JA,QA) ! tracer mixing ratio [kg/kg]

  ! diagnostic variables, defined at the cell center
  real(8), public, save :: VELZ(KA,IA,JA)    ! velocity w [m/s]
  real(8), public, save :: VELX(KA,IA,JA)    ! velocity u [m/s]
  real(8), public, save :: VELY(KA,IA,JA)    ! velocity v [m/s]
  real(8), public, save :: POTT(KA,IA,JA)    ! potential temperature [K]
  real(8), public, save :: QDRY(KA,IA,JA)    ! dry air mixig ratio [kg/kg]
  real(8), public, save :: PRES(KA,IA,JA)    ! pressure [Pa]
  real(8), public, save :: TEMP(KA,IA,JA)    ! temperature [K]
  real(8), public, save :: ENGT(KA,IA,JA)    ! total     energy [J/m3]
  real(8), public, save :: ENGP(KA,IA,JA)    ! potential energy [J/m3]
  real(8), public, save :: ENGK(KA,IA,JA)    ! kinetic   energy [J/m3]
  real(8), public, save :: ENGI(KA,IA,JA)    ! internal  energy [J/m3]

  character(len=FIO_HSHORT), public, save :: AP_NAME(5)
  character(len=FIO_HMID),   public, save :: AP_DESC(5)
  character(len=FIO_HSHORT), public, save :: AP_UNIT(5)

  data AP_NAME / 'DENS', &
                 'MOMZ', &
                 'MOMY', &
                 'MOMX', &
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

  character(len=IO_SYSCHR),  public, save :: ATMOS_TYPE_DYN    = 'NONE'
  character(len=IO_SYSCHR),  public, save :: ATMOS_TYPE_PHY_TB = 'NONE'
  character(len=IO_SYSCHR),  public, save :: ATMOS_TYPE_PHY_MP = 'NONE'
  character(len=IO_SYSCHR),  public, save :: ATMOS_TYPE_PHY_RD = 'NONE'

  logical,                   public, save :: ATMOS_sw_dyn
  logical,                   public, save :: ATMOS_sw_phy_tb
  logical,                   public, save :: ATMOS_sw_phy_mp
  logical,                   public, save :: ATMOS_sw_phy_rd
  logical,                   public, save :: ATMOS_sw_restart
  logical,                   public, save :: ATMOS_sw_check

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
  logical,                   private, save :: ATMOS_RESTART_IN_ALLOWMISSINGQ = .false.

  logical,                   private, save :: ATMOS_RESTART_CHECK            = .false.
  character(len=IO_FILECHR), private, save :: ATMOS_RESTART_CHECK_BASENAME   = 'restart_check'
  real(8),                   private, save :: ATMOS_RESTART_CHECK_CRITERION  = 1.D-6

  integer, private, save      :: ATMOS_HIST_sw (20)
  integer, private, save      :: ATMOS_PREP_sw (20)
  integer, private, save      :: ATMOS_MONIT_sw(20)
  integer, private, parameter :: I_VELZ =  1
  integer, private, parameter :: I_VELX =  2
  integer, private, parameter :: I_VELY =  3
  integer, private, parameter :: I_POTT =  4
  integer, private, parameter :: I_QDRY =  5
  integer, private, parameter :: I_RTOT =  6
  integer, private, parameter :: I_PRES =  7
  integer, private, parameter :: I_TEMP =  8
  integer, private, parameter :: I_RH   =  9
  integer, private, parameter :: I_QTOT = 10
  integer, private, parameter :: I_VOR  = 11
  integer, private, parameter :: I_DIV  = 12
  integer, private, parameter :: I_ENGP = 13
  integer, private, parameter :: I_ENGK = 14
  integer, private, parameter :: I_ENGI = 15
  integer, private, parameter :: I_ENGT = 16
  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Setup atmospheric variables
  !-----------------------------------------------------------------------------
  subroutine ATMOS_vars_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    use mod_const, only : &
       CPvap => CONST_CPvap, &
       CL    => CONST_CL,    &
       CI    => CONST_CI,    &
       CONST_UNDEF8
    use mod_history, only: &
       HIST_reg
    use mod_monitor, only: &
       MONIT_reg
    implicit none

    NAMELIST / PARAM_ATMOS / &
       ATMOS_TYPE_DYN,    &
       ATMOS_TYPE_PHY_TB, &
       ATMOS_TYPE_PHY_MP, &
       ATMOS_TYPE_PHY_RD

    NAMELIST / PARAM_ATMOS_VARS / &
       ATMOS_RESTART_IN_BASENAME,      &
       ATMOS_RESTART_IN_ALLOWMISSINGQ, &
       ATMOS_RESTART_OUTPUT,           &
       ATMOS_RESTART_OUT_BASENAME,     &
       ATMOS_RESTART_CHECK,            &
       ATMOS_RESTART_CHECK_BASENAME,   &
       ATMOS_RESTART_CHECK_CRITERION

    integer :: ierr
    integer :: iq
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

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** [ATMOS] selected components'

    if( IO_L ) write(IO_FID_LOG,*) 'Dynamics...'
    if ( ATMOS_TYPE_DYN == 'fent_fct' ) then
       if( IO_L ) write(IO_FID_LOG,*) '  Dynamical core   : Full-explicit No-terrain'
       if( IO_L ) write(IO_FID_LOG,*) '  Tracer advection : FCT limitter'
       ATMOS_sw_dyn = .true.
    elseif( ATMOS_TYPE_DYN == 'fent_pdfct' ) then
       if( IO_L ) write(IO_FID_LOG,*) '  Dynamical core   : Full-explicit No-terrain'
       if( IO_L ) write(IO_FID_LOG,*) '  Tracer advection : Positive-Difinite FCT limitter'
       ATMOS_sw_dyn = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '  Dynamical core   : NONE'
       if( IO_L ) write(IO_FID_LOG,*) '  Tracer advection : NONE'
       ATMOS_sw_dyn = .false.
    endif

    if( IO_L ) write(IO_FID_LOG,*) 'Physics...'
    if ( ATMOS_TYPE_PHY_TB == 'smagorinsky' ) then
       if( IO_L ) write(IO_FID_LOG,*) '  Sub-grid Turbulence : Smagorinsky'
       ATMOS_sw_phy_tb = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '  Sub-grid Turbulence : NONE'
       ATMOS_sw_phy_tb = .false.
    endif
    if ( ATMOS_TYPE_PHY_MP == 'NDW6' ) then
       if( IO_L ) write(IO_FID_LOG,*) '  Cloud Microphysics  : NDW6'
       ATMOS_sw_phy_mp = .true.
    elseif( ATMOS_TYPE_PHY_MP == 'kessler' ) then
       if( IO_L ) write(IO_FID_LOG,*) '  Cloud Microphysics  : Kessler'
       ATMOS_sw_phy_mp = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '  Cloud Microphysics  : NONE'
       ATMOS_sw_phy_mp = .false.
    endif
    if ( ATMOS_TYPE_PHY_RD == 'mstrnX' ) then
       if( IO_L ) write(IO_FID_LOG,*) '  Radiative transfer  : mstrnX'
       ATMOS_sw_phy_rd = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '  Radiative transfer  : NONE'
       ATMOS_sw_phy_rd = .false.
    endif

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
    if( IO_L ) write(IO_FID_LOG,*) '***       | VARNAME|DESCRIPTION', &
    '                                                     [UNIT            ]'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,i3,A,A8,5(A))') &
    '*** NO.',1,'|',trim(AP_NAME(1)),'|', AP_DESC(1),'[', AP_UNIT(1),']'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,i3,A,A8,5(A))') &
    '*** NO.',2,'|',trim(AP_NAME(2)),'|', AP_DESC(2),'[', AP_UNIT(2),']'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,i3,A,A8,5(A))') &
    '*** NO.',3,'|',trim(AP_NAME(3)),'|', AP_DESC(3),'[', AP_UNIT(3),']'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,i3,A,A8,5(A))') &
    '*** NO.',4,'|',trim(AP_NAME(4)),'|', AP_DESC(4),'[', AP_UNIT(4),']'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,i3,A,A8,5(A))') &
    '*** NO.',5,'|',trim(AP_NAME(5)),'|', AP_DESC(5),'[', AP_UNIT(5),']'
    do iq = 1, QA
       if( IO_L ) write(IO_FID_LOG,'(1x,A,i3,A,A8,5(A))')  &
       '*** NO.',5+iq,'|',trim(AQ_NAME(iq)),'|', AQ_DESC(iq),'[', AQ_UNIT(iq),']'
    enddo

    if( IO_L ) write(IO_FID_LOG,*) 'Output...'
    if ( ATMOS_RESTART_OUTPUT ) then
       if( IO_L ) write(IO_FID_LOG,*) '  Restart output : YES'
       ATMOS_sw_restart = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '  Restart output : NO'
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

    ATMOS_HIST_sw (:) = -1
    ATMOS_MONIT_sw(:) = -1
    ATMOS_PREP_sw (:) = -1

    call HIST_reg( ATMOS_HIST_sw(I_VELZ), 'W',    'velocity w',            'm/s',   '3D' )
    call HIST_reg( ATMOS_HIST_sw(I_VELX), 'U',    'velocity u',            'm/s',   '3D' )
    call HIST_reg( ATMOS_HIST_sw(I_VELY), 'V',    'velocity v',            'm/s',   '3D' )
    call HIST_reg( ATMOS_HIST_sw(I_POTT), 'PT',   'potential temp.',       'K',     '3D' )
    call HIST_reg( ATMOS_HIST_sw(I_QDRY), 'QDRY', 'dry air',               'kg/kg', '3D' )
    call HIST_reg( ATMOS_HIST_sw(I_RTOT), 'RTOT', 'Total gas constant',    'kg/kg', '3D' )
    call HIST_reg( ATMOS_HIST_sw(I_PRES), 'PRES', 'pressure',              'Pa',    '3D' )
    call HIST_reg( ATMOS_HIST_sw(I_TEMP), 'T',    'temperature',           'K',     '3D' )
    call HIST_reg( ATMOS_HIST_sw(I_RH  ), 'RH',   'relative humidity',     '%',     '3D' )
    call HIST_reg( ATMOS_HIST_sw(I_QTOT), 'QTOT', 'total hydrometeors',    'kg/kg', '3D' )
    call HIST_reg( ATMOS_HIST_sw(I_VOR ), 'VOR',  'vertical vorticity',    '1/s',   '3D' )
    call HIST_reg( ATMOS_HIST_sw(I_DIV ), 'DIV',  'horizontal divergence', '1/s',   '3D' )

    call HIST_reg( ATMOS_HIST_sw(I_ENGP), 'ENGP', 'potential energy',      'J/m3',  '3D' )
    call HIST_reg( ATMOS_HIST_sw(I_ENGK), 'ENGK', 'kinetic energy',        'J/m3',  '3D' )
    call HIST_reg( ATMOS_HIST_sw(I_ENGI), 'ENGI', 'internal energy',       'J/m3',  '3D' )
    call MONIT_reg( ATMOS_MONIT_sw(I_ENGP), 'ENGP', 'potential energy', 'J', '3D' )
    call MONIT_reg( ATMOS_MONIT_sw(I_ENGK), 'ENGK', 'kinetic   energy', 'J', '3D' )
    call MONIT_reg( ATMOS_MONIT_sw(I_ENGI), 'ENGI', 'internal  energy', 'J', '3D' )
    call MONIT_reg( ATMOS_MONIT_sw(I_ENGT), 'ENGT', 'total     energy', 'J', '3D' )

    if ( ATMOS_HIST_sw(I_VELZ) > 0 ) then
       ATMOS_PREP_sw(I_VELZ) = 1
    endif
    if ( ATMOS_HIST_sw(I_VELX) > 0 ) then
       ATMOS_PREP_sw(I_VELX) = 1
    endif
    if ( ATMOS_HIST_sw(I_VELY) > 0 ) then
       ATMOS_PREP_sw(I_VELY) = 1
    endif
    if ( ATMOS_HIST_sw(I_POTT) > 0 ) then
       ATMOS_PREP_sw(I_POTT) = 1
    endif
    if ( ATMOS_HIST_sw(I_QDRY) > 0 ) then
       ATMOS_PREP_sw(I_QDRY) = 1
    endif
    if ( ATMOS_HIST_sw(I_RTOT) > 0 ) then
       ATMOS_PREP_sw(I_QDRY) = 1
       ATMOS_PREP_sw(I_RTOT) = 1
    endif
    if ( ATMOS_HIST_sw(I_PRES) > 0 ) then
       ATMOS_PREP_sw(I_QDRY) = 1
       ATMOS_PREP_sw(I_RTOT) = 1
       ATMOS_PREP_sw(I_PRES) = 1
    endif
    if ( ATMOS_HIST_sw(I_TEMP) > 0 ) then
       ATMOS_PREP_sw(I_QDRY) = 1
       ATMOS_PREP_sw(I_RTOT) = 1
       ATMOS_PREP_sw(I_PRES) = 1
       ATMOS_PREP_sw(I_TEMP) = 1
    endif
    if ( ATMOS_HIST_sw(I_RH) > 0 ) then
       ATMOS_PREP_sw(I_QDRY) = 1
       ATMOS_PREP_sw(I_RTOT) = 1
       ATMOS_PREP_sw(I_PRES) = 1
       ATMOS_PREP_sw(I_TEMP) = 1
       ATMOS_PREP_sw(I_RH  ) = 1
    endif
    if ( ATMOS_HIST_sw (I_QTOT) > 0 ) then
       ATMOS_PREP_sw(I_QTOT) = 1
    endif
    if ( ATMOS_HIST_sw (I_VOR) > 0 ) then
       ATMOS_PREP_sw(I_VOR) = 1
    endif

    if (      ATMOS_HIST_sw (I_ENGP) > 0 &
         .OR. ATMOS_MONIT_sw(I_ENGP) > 0 ) then
       ATMOS_PREP_sw(I_ENGP) = 1
    endif
    if (      ATMOS_HIST_sw (I_ENGK) > 0 &
         .OR. ATMOS_MONIT_sw(I_ENGK) > 0 ) then
       ATMOS_PREP_sw(I_VELZ) = 1
       ATMOS_PREP_sw(I_VELX) = 1
       ATMOS_PREP_sw(I_VELY) = 1
       ATMOS_PREP_sw(I_ENGK) = 1
    endif
    if (      ATMOS_HIST_sw (I_ENGI) > 0 &
         .OR. ATMOS_MONIT_sw(I_ENGI) > 0 ) then
       ATMOS_PREP_sw(I_QDRY) = 1
       ATMOS_PREP_sw(I_RTOT) = 1
       ATMOS_PREP_sw(I_PRES) = 1
       ATMOS_PREP_sw(I_TEMP) = 1
       ATMOS_PREP_sw(I_ENGI) = 1
    endif
    if (      ATMOS_HIST_sw (I_ENGT) > 0 &
         .OR. ATMOS_MONIT_sw(I_ENGT) > 0 ) then
       ATMOS_PREP_sw(I_ENGP) = 1
       ATMOS_PREP_sw(I_ENGK) = 1
       ATMOS_PREP_sw(I_ENGI) = 1
       ATMOS_PREP_sw(I_ENGT) = 1
    endif

    return
  end subroutine ATMOS_vars_setup

  !-----------------------------------------------------------------------------
  !> fill HALO region of atmospheric variables
  !-----------------------------------------------------------------------------
  subroutine ATMOS_vars_fillhalo
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    integer :: i, j, iq
    !---------------------------------------------------------------------------

    ! fill KHALO
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
  !-----------------------------------------------------------------------------
  subroutine ATMOS_vars_restart_read
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait
    use mod_fileio, only: &
       FIO_input
    implicit none

    real(8) :: restart_atmos(KMAX,IMAX,JMAX) !> restart file (no HALO)

    character(len=IO_FILECHR) :: bname
    character(len=8)          :: lname

    integer :: iq
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input restart file (atmos) ***'

    bname = ATMOS_RESTART_IN_BASENAME
    write(lname,'(A,I4.4)') 'ZDEF', KMAX

    call FIO_input( restart_atmos(:,:,:), bname, 'DENS', lname, 1, KMAX, 1 )
    DENS(KS:KE,IS:IE,JS:JE) = restart_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FIO_input( restart_atmos(:,:,:), bname, 'MOMZ', lname, 1, KMAX, 1 )
    MOMZ(KS:KE,IS:IE,JS:JE) = restart_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FIO_input( restart_atmos(:,:,:), bname, 'MOMX', lname, 1, KMAX, 1 )
    MOMX(KS:KE,IS:IE,JS:JE) = restart_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FIO_input( restart_atmos(:,:,:), bname, 'MOMY', lname, 1, KMAX, 1 )
    MOMY(KS:KE,IS:IE,JS:JE) = restart_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FIO_input( restart_atmos(:,:,:), bname, 'RHOT', lname, 1, KMAX, 1 )
    RHOT(KS:KE,IS:IE,JS:JE) = restart_atmos(1:KMAX,1:IMAX,1:JMAX)

    do iq = 1, QA
       call FIO_input( restart_atmos(:,:,:), bname, AQ_NAME(iq), lname, 1, KMAX, 1 )
       QTRC(KS:KE,IS:IE,JS:JE,iq) = restart_atmos(1:KMAX,1:IMAX,1:JMAX)
    enddo

    call ATMOS_vars_fillhalo

    call ATMOS_vars_total

    return
  end subroutine ATMOS_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Write restart of atmospheric variables
  !-----------------------------------------------------------------------------
  subroutine ATMOS_vars_restart_write
    use mod_time, only: &
       NOWSEC => TIME_NOWSEC
    use mod_fileio_h, only: &
       FIO_REAL8
    use mod_fileio, only: &
       FIO_output
    implicit none

    real(8) :: restart_atmos(KMAX,IMAX,JMAX) !> restart file (no HALO)

    character(len=IO_FILECHR) :: bname
    character(len=FIO_HMID)   :: desc
    character(len=8)          :: lname

    integer :: iq
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Output restart file (atmos) ***'

    write(bname,'(A,A,F15.3)') trim(ATMOS_RESTART_OUT_BASENAME), '_', NOWSEC
    desc  = 'SCALE3 PROGNOSTIC VARS.'
    write(lname,'(A,I4.4)') 'ZDEF', KMAX

    restart_atmos(1:KMAX,1:IMAX,1:JMAX) = DENS(KS:KE,IS:IE,JS:JE)
    call FIO_output( restart_atmos(:,:,:), bname, desc, '',       &
                     'DENS', AP_DESC(1), '', AP_UNIT(1),          &
                     FIO_REAL8, lname, 1, KMAX, 1, NOWSEC, NOWSEC )

    restart_atmos(1:KMAX,1:IMAX,1:JMAX) = MOMZ(KS:KE,IS:IE,JS:JE)
    call FIO_output( restart_atmos(:,:,:), bname, desc, '',       &
                     'MOMZ', AP_DESC(2), '', AP_UNIT(2),          &
                     FIO_REAL8, lname, 1, KMAX, 1, NOWSEC, NOWSEC )

    restart_atmos(1:KMAX,1:IMAX,1:JMAX) = MOMX(KS:KE,IS:IE,JS:JE)
    call FIO_output( restart_atmos(:,:,:), bname, desc, '',       &
                     'MOMX', AP_DESC(3), '', AP_UNIT(3),          &
                     FIO_REAL8, lname, 1, KMAX, 1, NOWSEC, NOWSEC )

    restart_atmos(1:KMAX,1:IMAX,1:JMAX) = MOMY(KS:KE,IS:IE,JS:JE)
    call FIO_output( restart_atmos(:,:,:), bname, desc, '',       &
                     'MOMY', AP_DESC(4), '', AP_UNIT(4),          &
                     FIO_REAL8, lname, 1, KMAX, 1, NOWSEC, NOWSEC )

    restart_atmos(1:KMAX,1:IMAX,1:JMAX) = RHOT(KS:KE,IS:IE,JS:JE)
    call FIO_output( restart_atmos(:,:,:), bname, desc, '',       &
                     'RHOT', AP_DESC(5), '', AP_UNIT(5),          &
                     FIO_REAL8, lname, 1, KMAX, 1, NOWSEC, NOWSEC )

    do iq = 1, QA
       restart_atmos(1:KMAX,1:IMAX,1:JMAX) = QTRC(KS:KE,IS:IE,JS:JE,iq)

       call FIO_output( restart_atmos(:,:,:), bname, desc, '',       &
                        AQ_NAME(iq), AQ_DESC(iq), '', AQ_UNIT(iq),   &
                        FIO_REAL8, lname, 1, KMAX, 1, NOWSEC, NOWSEC )
    enddo

    call ATMOS_vars_total

    return
  end subroutine ATMOS_vars_restart_write

  !-----------------------------------------------------------------------------
  !> Check and compare between last data and sample data
  !-----------------------------------------------------------------------------
  subroutine ATMOS_vars_restart_check
    use mod_process, only: &
       PRC_myrank
    use mod_fileio, only: &
       FIO_input
    implicit none

    real(8) :: DENS_check(KA,IA,JA)    ! Density    [kg/m3]
    real(8) :: MOMZ_check(KA,IA,JA)    ! momentum z [kg/s/m2]
    real(8) :: MOMX_check(KA,IA,JA)    ! momentum x [kg/s/m2]
    real(8) :: MOMY_check(KA,IA,JA)    ! momentum y [kg/s/m2]
    real(8) :: RHOT_check(KA,IA,JA)    ! DENS * POTT [K*kg/m3]
    real(8) :: QTRC_check(KA,IA,JA,QA) ! tracer mixing ratio [kg/kg]

    real(8) :: restart_atmos(KMAX,IMAX,JMAX) !> restart file (no HALO)

    character(len=IO_FILECHR) :: bname
    character(len=8)          :: lname

    logical :: datacheck
    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    write(*,*) 'Compare last Data with ', trim(ATMOS_RESTART_CHECK_BASENAME), 'on PE=', PRC_myrank
    write(*,*) '*** criterion = ', ATMOS_RESTART_CHECK_CRITERION
    datacheck = .true.

    bname = ATMOS_RESTART_CHECK_BASENAME
    write(lname,'(A,I4.4)') 'ZDEF', KMAX

    call FIO_input( restart_atmos(:,:,:), bname, 'DENS', lname, 1, KMAX, 1 )
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

    call FIO_input( restart_atmos(:,:,:), bname, 'MOMZ', lname, 1, KMAX, 1 )
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

    call FIO_input( restart_atmos(:,:,:), bname, 'MOMX', lname, 1, KMAX, 1 )
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

    call FIO_input( restart_atmos(:,:,:), bname, 'MOMY', lname, 1, KMAX, 1 )
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

    call FIO_input( restart_atmos(:,:,:), bname, 'RHOT', lname, 1, KMAX, 1 )
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
       call FIO_input( restart_atmos(:,:,:), bname, AQ_NAME(iq), lname, 1, KMAX, 1 )
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

    return
  end subroutine ATMOS_vars_restart_check

  !-----------------------------------------------------------------------------
  !> History output set for prognostic variables
  !-----------------------------------------------------------------------------
  subroutine ATMOS_vars_history
    use mod_const, only : &
       GRAV   => CONST_GRAV,   &
       Rdry   => CONST_Rdry,   &
       Rvap   => CONST_Rvap,   &
       CVdry  => CONST_CVdry,  &
       CPovCV => CONST_CPovCV, &
       P00    => CONST_PRE00
    use mod_time, only: &
       TIME_DTSEC
    use mod_grid, only : &
       CZ   => GRID_CZ,   &
       RCDX => GRID_RCDX, &
       RCDY => GRID_RCDY
    use mod_history, only: &
       HIST_put, &
       HIST_in
    use mod_monitor, only: &
       MONIT_put, &
       MONIT_in
    use mod_atmos_thermodyn, only: &
       CVw => AQ_CV
    use mod_atmos_saturation, only: &
       saturation_qsat_water => ATMOS_SATURATION_qsat_water
    implicit none

    real(8) :: RTOT (KA,IA,JA)
    real(8) :: QSAT (KA,IA,JA)
    real(8) :: RH   (KA,IA,JA)
    real(8) :: QTOT (KA,IA,JA)
    real(8) :: VELXH(KA,IA,JA)
    real(8) :: VELYH(KA,IA,JA)
    real(8) :: VOR  (KA,IA,JA)

    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    call HIST_in( DENS(:,:,:), 'DENS', AP_DESC(1), AP_UNIT(1), '3D', TIME_DTSEC )
    call HIST_in( MOMZ(:,:,:), 'MOMZ', AP_DESC(2), AP_UNIT(2), '3D', TIME_DTSEC )
    call HIST_in( MOMX(:,:,:), 'MOMX', AP_DESC(3), AP_UNIT(3), '3D', TIME_DTSEC )
    call HIST_in( MOMY(:,:,:), 'MOMY', AP_DESC(4), AP_UNIT(4), '3D', TIME_DTSEC )
    call HIST_in( RHOT(:,:,:), 'RHOT', AP_DESC(5), AP_UNIT(5), '3D', TIME_DTSEC )
    do iq = 1, QA
       call HIST_in( QTRC(:,:,:,iq), AQ_NAME(iq), AQ_DESC(iq), AQ_UNIT(iq), '3D', TIME_DTSEC )
    enddo

    call MONIT_in( DENS(:,:,:), 'DENS', AP_DESC(1), AP_UNIT(1), '3D' )
    call MONIT_in( MOMZ(:,:,:), 'MOMZ', AP_DESC(2), AP_UNIT(2), '3D' )
    call MONIT_in( MOMX(:,:,:), 'MOMX', AP_DESC(3), AP_UNIT(3), '3D' )
    call MONIT_in( MOMY(:,:,:), 'MOMY', AP_DESC(4), AP_UNIT(4), '3D' )
    call MONIT_in( RHOT(:,:,:), 'RHOT', AP_DESC(5), AP_UNIT(5), '3D' )
    do iq = 1, QA
       call MONIT_in( QTRC(:,:,:,iq), AQ_NAME(iq), AQ_DESC(iq), AQ_UNIT(iq), '3D' )
    enddo

    if ( ATMOS_PREP_sw(I_VELZ) > 0 ) then
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          VELZ(k,i,j) = 0.5D0 * ( MOMZ(k-1,i,j)+MOMZ(k,i,j) ) / DENS(k,i,j)
       enddo
       enddo
       enddo
    endif

    if ( ATMOS_PREP_sw(I_VELX) > 0 ) then
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          VELX(k,i,j) = 0.5D0 * ( MOMX(k,i-1,j)+MOMX(k,i,j) ) / DENS(k,i,j)
       enddo
       enddo
       enddo
    endif

    if ( ATMOS_PREP_sw(I_VELY) > 0 ) then
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          VELY(k,i,j) = 0.5D0 * ( MOMY(k,i,j-1)+MOMY(k,i,j) ) / DENS(k,i,j)
       enddo
       enddo
       enddo
    endif

    if ( ATMOS_PREP_sw(I_POTT) > 0 ) then
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          POTT(k,i,j) = RHOT(k,i,j) / DENS(k,i,j)
       enddo
       enddo
       enddo
    endif

    if ( ATMOS_PREP_sw(I_QDRY) > 0 ) then
       do j  = JS, JE
       do i  = IS, IE
       do k  = KS, KE
          QDRY(k,i,j) = 1.D0
       enddo
       enddo
       enddo

       do iq = QQS, QQE
       do j  = JS, JE
       do i  = IS, IE
       do k  = KS, KE
          QDRY(k,i,j) = QDRY(k,i,j) - QTRC(k,i,j,iq)
       enddo
       enddo
       enddo
       enddo
    endif

    if ( ATMOS_PREP_sw(I_RTOT) > 0 ) then
       do j  = JS, JE
       do i  = IS, IE
       do k  = KS, KE
          RTOT(k,i,j) = Rdry*QDRY(k,i,j) + Rvap*QTRC(k,i,j,I_QV)
       enddo
       enddo
       enddo
    endif

    if ( ATMOS_PREP_sw(I_PRES) > 0 ) then
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          PRES(k,i,j) = P00 * ( RHOT(k,i,j) * RTOT(k,i,j) / P00 )**CPovCV
       enddo
       enddo
       enddo
    endif

    if ( ATMOS_PREP_sw(I_TEMP) > 0 ) then
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          TEMP(k,i,j) = PRES(k,i,j) / ( DENS(k,i,j) * RTOT(k,i,j) )
       enddo
       enddo
       enddo
    endif

    if ( ATMOS_HIST_sw(I_RH) > 0 ) then
       call saturation_qsat_water( QSAT(:,:,:), TEMP(:,:,:), PRES(:,:,:) )

       do j  = JS, JE
       do i  = IS, IE
       do k  = KS, KE
          RH(k,i,j) = QTRC(k,i,j,I_QV) / QSAT(k,i,j) * 1.D2
       enddo
       enddo
       enddo
    endif

    if ( ATMOS_PREP_sw(I_QTOT) > 0 ) then
       do j  = JS, JE
       do i  = IS, IE
       do k  = KS, KE
          QTOT(k,i,j) = 0.D0
       enddo
       enddo
       enddo

       if ( iq >= 2 ) then
          do iq = QQS+1, QQE
          do j  = JS, JE
          do i  = IS, IE
          do k  = KS, KE
             QTOT(k,i,j) = QTOT(k,i,j) + QTRC(k,i,j,iq)
          enddo
          enddo
          enddo
          enddo
       endif
    endif

    if ( ATMOS_PREP_sw(I_VOR) > 0 ) then
       ! at u, v, layer
       do j = JS-1, JE
       do i = IS-1, IE
       do k = KS, KE
          VELXH(k,i,j) = 2.D0 * ( MOMX(k,i,j)+MOMX(k,i,j+1) )                             &
                              / ( DENS(k,i,j)+DENS(k,i+1,j)+DENS(k,i,j+1)+DENS(k,i+1,j+1) )
       enddo
       enddo
       enddo

       ! at u, v, layer
       do j = JS-1, JE
       do i = IS-1, IE
       do k = KS, KE
          VELYH(k,i,j) = 2.D0 * ( MOMY(k,i,j)+MOMY(k,i+1,j) )                             &
                              / ( DENS(k,i,j)+DENS(k,i+1,j)+DENS(k,i,j+1)+DENS(k,i+1,j+1) )
       enddo
       enddo
       enddo

       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          VOR(k,i,j) = 0.5D0 * ( ( VELYH(k,i,j  )-VELYH(k,i-1,j  ) ) * RCDX(i) &
                               + ( VELYH(k,i,j-1)-VELYH(k,i-1,j-1) ) * RCDX(i) &
                               - ( VELXH(k,i  ,j)-VELXH(k,i  ,j-1) ) * RCDY(j) &
                               - ( VELXH(k,i-1,j)-VELXH(k,i-1,j-1) ) * RCDY(j) )
       enddo
       enddo
       enddo
    endif

    if ( ATMOS_PREP_sw(I_ENGP) > 0 ) then
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          ENGP(k,i,j) = DENS(k,i,j) * GRAV * CZ(k)
       enddo
       enddo
       enddo
    endif

    if ( ATMOS_PREP_sw(I_ENGK) > 0 ) then
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          ENGK(k,i,j) = 0.5D0 * DENS(k,i,j)    &
                              * VELZ(k,i,j)**2 &
                              * VELX(k,i,j)**2 &
                              * VELY(k,i,j)**2
       enddo
       enddo
       enddo
    endif

    if ( ATMOS_PREP_sw(I_ENGI) > 0 ) then
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          ENGI(k,i,j) = DENS(k,i,j) * QDRY(k,i,j) * TEMP(k,i,j) * CVdry
       enddo
       enddo
       enddo

       do iq = QQS, QQE
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

    if ( ATMOS_HIST_sw(I_ENGT) > 0 ) then
       do j  = JS, JE
       do i  = IS, IE
       do k  = KS, KE
          ENGT(k,i,j) = ENGP(k,i,j) + ENGK(k,i,j) + ENGI(k,i,j)
       enddo
       enddo
       enddo
    endif

    if( ATMOS_HIST_sw(I_VELZ) > 0 ) call HIST_put( ATMOS_HIST_sw(I_VELZ), VELZ(:,:,:), TIME_DTSEC )
    if( ATMOS_HIST_sw(I_VELX) > 0 ) call HIST_put( ATMOS_HIST_sw(I_VELX), VELX(:,:,:), TIME_DTSEC )
    if( ATMOS_HIST_sw(I_VELY) > 0 ) call HIST_put( ATMOS_HIST_sw(I_VELY), VELY(:,:,:), TIME_DTSEC )
    if( ATMOS_HIST_sw(I_POTT) > 0 ) call HIST_put( ATMOS_HIST_sw(I_POTT), POTT(:,:,:), TIME_DTSEC )
    if( ATMOS_HIST_sw(I_QDRY) > 0 ) call HIST_put( ATMOS_HIST_sw(I_QDRY), QDRY(:,:,:), TIME_DTSEC )
    if( ATMOS_HIST_sw(I_PRES) > 0 ) call HIST_put( ATMOS_HIST_sw(I_PRES), PRES(:,:,:), TIME_DTSEC )
    if( ATMOS_HIST_sw(I_TEMP) > 0 ) call HIST_put( ATMOS_HIST_sw(I_TEMP), TEMP(:,:,:), TIME_DTSEC )
    if( ATMOS_HIST_sw(I_RH  ) > 0 ) call HIST_put( ATMOS_HIST_sw(I_RH  ), RH  (:,:,:), TIME_DTSEC )
    if( ATMOS_HIST_sw(I_QTOT) > 0 ) call HIST_put( ATMOS_HIST_sw(I_QTOT), QTOT(:,:,:), TIME_DTSEC )
    if( ATMOS_HIST_sw(I_VOR ) > 0 ) call HIST_put( ATMOS_HIST_sw(I_VOR),  VOR (:,:,:), TIME_DTSEC )

    if( ATMOS_HIST_sw(I_ENGP) > 0 ) call HIST_put( ATMOS_HIST_sw(I_ENGP), ENGP(:,:,:), TIME_DTSEC )
    if( ATMOS_HIST_sw(I_ENGK) > 0 ) call HIST_put( ATMOS_HIST_sw(I_ENGK), ENGK(:,:,:), TIME_DTSEC )
    if( ATMOS_HIST_sw(I_ENGI) > 0 ) call HIST_put( ATMOS_HIST_sw(I_ENGI), ENGI(:,:,:), TIME_DTSEC )

    if( ATMOS_MONIT_sw(I_ENGP) > 0 ) call MONIT_put( ATMOS_MONIT_sw(I_ENGP), ENGP(:,:,:) )
    if( ATMOS_MONIT_sw(I_ENGK) > 0 ) call MONIT_put( ATMOS_MONIT_sw(I_ENGK), ENGK(:,:,:) )
    if( ATMOS_MONIT_sw(I_ENGI) > 0 ) call MONIT_put( ATMOS_MONIT_sw(I_ENGI), ENGI(:,:,:) )
    if( ATMOS_MONIT_sw(I_ENGT) > 0 ) call MONIT_put( ATMOS_MONIT_sw(I_ENGT), ENGT(:,:,:) )

    return
  end subroutine ATMOS_vars_history

  !-----------------------------------------------------------------------------
  !> Budget monitor of atmosphere
  !-----------------------------------------------------------------------------
  subroutine ATMOS_vars_total
    use mod_const, only : &
       GRAV   => CONST_GRAV,   &
       Rdry   => CONST_Rdry,   &
       CVdry  => CONST_CVdry,  &
       CPovCV => CONST_CPovCV, &
       P00    => CONST_PRE00
    use mod_grid, only : &
       CZ   => GRID_CZ
    use mod_comm, only: &
       COMM_total_doreport, &
       COMM_total
    use mod_atmos_thermodyn, only: &
       CVw => AQ_CV
    implicit none

    integer :: i, j, k, iq
    !---------------------------------------------------------------------------

    if ( COMM_total_doreport ) then

       call COMM_total( DENS(:,:,:), AP_NAME(1) )
       call COMM_total( MOMZ(:,:,:), AP_NAME(2) )
       call COMM_total( MOMX(:,:,:), AP_NAME(3) )
       call COMM_total( MOMY(:,:,:), AP_NAME(4) )
       call COMM_total( RHOT(:,:,:), AP_NAME(5) )
       do iq = 1, QA
          call COMM_total( QTRC(:,:,:,iq), AQ_NAME(iq) )
       enddo

       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          VELZ(k,i,j) = 0.5D0 * ( MOMZ(k-1,i,j)+MOMZ(k,i,j) ) / DENS(k,i,j)
          VELX(k,i,j) = 0.5D0 * ( MOMX(k,i-1,j)+MOMX(k,i,j) ) / DENS(k,i,j)
          VELY(k,i,j) = 0.5D0 * ( MOMY(k,i,j-1)+MOMY(k,i,j) ) / DENS(k,i,j)
          PRES(k,i,j) = P00 * ( RHOT(k,i,j) * Rdry / P00 )**CPovCV
          TEMP(k,i,j) = PRES(k,i,j) / ( DENS(k,i,j) * Rdry )
          ENGP(k,i,j) = DENS(k,i,j) * GRAV * CZ(k)
          ENGK(k,i,j) = 0.5D0 * DENS(k,i,j)    &
                              * VELZ(k,i,j)**2 &
                              * VELX(k,i,j)**2 &
                              * VELY(k,i,j)**2
          QDRY(k,i,j) = 1.D0
          do iq = QQS, QQE
             QDRY(k,i,j) = QDRY(k,i,j) - QTRC(k,i,j,iq)
          enddo
          ENGI(k,i,j) = DENS(k,i,j) * QDRY(k,i,j) * TEMP(k,i,j) * CVdry
          do iq = QQS, QQE
             ENGI(k,i,j) = ENGI(k,i,j) &
                         + DENS(k,i,j) * QTRC(k,i,j,iq) * TEMP(k,i,j) * CVw(iq)
          enddo
          ENGT(k,i,j) = ENGP(k,i,j) + ENGK(k,i,j) + ENGI(k,i,j)
       enddo
       enddo
       enddo

       call COMM_total( ENGT(:,:,:), 'ENGT    ' )
       call COMM_total( ENGP(:,:,:), 'ENGP    ' )
       call COMM_total( ENGK(:,:,:), 'ENGK    ' )
       call COMM_total( ENGI(:,:,:), 'ENGI    ' )

    endif

    return
  end subroutine ATMOS_vars_total

end module mod_atmos_vars
