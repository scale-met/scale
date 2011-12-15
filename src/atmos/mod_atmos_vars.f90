!-------------------------------------------------------------------------------
!> module Atmospheric Variables
!!
!! @par Description
!!          Container for atmospheric variables
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2011-11-11 (H.Yashiro) [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_vars
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_stdio, only: &
     IO_SYSCHR, &
     IO_FILECHR
  use mod_fileio_h, only: &
     FIO_HSHORT, &
     FIO_HMID,   &
     FIO_REAL8
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_vars_setup
  public :: ATMOS_vars_restart_read
  public :: ATMOS_vars_restart_write
  public :: ATMOS_vars_restart_check
  public :: ATMOS_vars_put
  public :: ATMOS_vars_get
  public :: ATMOS_vars_getall
  public :: ATMOS_DMP2PVT
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer,                   public,              save :: A_VA      ! Number of Tracers + 5
  integer,                   public,              save :: A_QA      ! Number of Tracers
  character(len=FIO_HSHORT), public, allocatable, save :: A_NAME(:)
  character(len=FIO_HMID),   public, allocatable, save :: A_DESC(:)
  character(len=FIO_HSHORT), public, allocatable, save :: A_UNIT(:)

  integer, public, parameter :: A_QWA_MAX = 6
  real(8), public,      save :: A_CPw(A_QWA_MAX) = 0.D0 ! CP for each hydrometeors
  real(8), public,      save :: A_CVw(A_QWA_MAX) = 0.D0 ! CV for each hydrometeors

  integer, public,      save :: A_QWA =  0 ! Total number of all hydrometeors
  integer, public,      save :: A_QWS = -1 ! Start index for all hydrometeors
  integer, public,      save :: A_QWE = -1 ! End   index for all hydrometeors
  integer, public,      save :: A_QLS = -1 ! Start index for liquids
  integer, public,      save :: A_QLE = -1 ! End   index for liquids
  integer, public,      save :: A_QSS = -1 ! Start index for solids
  integer, public,      save :: A_QSE = -1 ! End   index for solids
  integer, public,      save :: A_NWS = -1 ! Start index for number of liquids, solids
  integer, public,      save :: A_NWE = -1 ! End   index for number of liquids, solids

  integer, public,      save :: I_QV  = -1
  integer, public,      save :: I_QC  = -1
  integer, public,      save :: I_QR  = -1
  integer, public,      save :: I_QI  = -1
  integer, public,      save :: I_QS  = -1
  integer, public,      save :: I_QG  = -1
  integer, public,      save :: I_NC  = -1
  integer, public,      save :: I_NR  = -1
  integer, public,      save :: I_NI  = -1
  integer, public,      save :: I_NS  = -1
  integer, public,      save :: I_NG  = -1

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

  real(8), public, allocatable, save :: atmos_var(:,:,:,:)      !> prognostics container (with HALO)
  real(8), public, allocatable, save :: atmos_diagvar(:,:,:,:)  !> diagnostics container (with HALO)

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

  logical,                   private, save :: ATMOS_RESTART_CHECK            = .true.
  character(len=IO_FILECHR), private, save :: ATMOS_RESTART_CHECK_BASENAME   = 'restart_check'
  real(8),                   private, save :: ATMOS_RESTART_CHECK_CRITERION  = 1.D-6

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Setup atmospheric variables
  !-----------------------------------------------------------------------------
  subroutine ATMOS_vars_setup
    use mod_stdio, only: &
       IO_FID_CONF, &
       IO_FID_LOG,  &
       IO_L
    use mod_process, only: &
       PRC_MPIstop
    use mod_const, only : &
       CPvap => CONST_CPvap, &
       CVvap => CONST_CVvap, &
       CL    => CONST_CL,    &
       CI    => CONST_CI,    &
       CONST_UNDEF8
    use mod_grid, only: &
       IA => GRID_IA, &
       JA => GRID_JA, &
       KA => GRID_KA
    implicit none

    NAMELIST / PARAM_ATMOS / &
       ATMOS_TYPE_DYN,    &
       ATMOS_TYPE_PHY_TB, &
       ATMOS_TYPE_PHY_MP, &
       ATMOS_TYPE_PHY_RD

    integer                   :: ATMOS_QTRC_NMAX    = 0
    character(len=FIO_HSHORT) :: ATMOS_QTRC_VARNAME = ''
    character(len=IO_FILECHR) :: ATMOS_QTRC_VARDESC = ''
    character(len=IO_FILECHR) :: ATMOS_QTRC_VARUNIT = ''

    NAMELIST / PARAM_ATMOS_VARS / &
       ATMOS_QTRC_NMAX,                &
       ATMOS_QTRC_VARNAME,             &
       ATMOS_QTRC_VARDESC,             &
       ATMOS_QTRC_VARUNIT,             &
       ATMOS_RESTART_OUTPUT,           &
       ATMOS_RESTART_IN_BASENAME,      &
       ATMOS_RESTART_OUT_BASENAME,     &
       ATMOS_RESTART_IN_ALLOWMISSINGQ, &
       ATMOS_RESTART_CHECK,            &
       ATMOS_RESTART_CHECK_BASENAME,   &
       ATMOS_RESTART_CHECK_CRITERION

    integer :: ierr
    integer :: iv
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

    if( IO_L ) write(IO_FID_LOG,*) '*** Dynamics...'
    if ( ATMOS_TYPE_DYN == 'fent_fct' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Dynamical core   : Full-explicit No-terrain'
       if( IO_L ) write(IO_FID_LOG,*) '*** Tracer advection : FCT limitter'
       ATMOS_sw_dyn = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Dynamical core   : NONE'
       if( IO_L ) write(IO_FID_LOG,*) '*** Tracer advection : NONE'
       ATMOS_sw_dyn = .false.
    endif

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics...'
    if ( ATMOS_TYPE_PHY_TB == 'smagorinsky' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Sub-grid Turbulence : Smagorinsky'
       ATMOS_sw_phy_tb = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Sub-grid Turbulence : NONE'
       ATMOS_sw_phy_tb = .false.
    endif
    if ( ATMOS_TYPE_PHY_MP == 'NDW6' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Cloud Microphysics  : NDW6'
       ATMOS_sw_phy_mp = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Cloud Microphysics  : NONE'
       ATMOS_sw_phy_mp = .false.
    endif
    if ( ATMOS_TYPE_PHY_RD == 'mstrnX' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Radiative transfer  : mstrnX'
       ATMOS_sw_phy_rd = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Radiative transfer  : NONE'
       ATMOS_sw_phy_rd = .false.
    endif

    if( IO_L ) write(IO_FID_LOG,*) '*** Output...'
    if ( ATMOS_RESTART_OUTPUT ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart output : YES'
       ATMOS_sw_restart = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart output : NO'
       ATMOS_sw_restart = .false.
    endif
    if ( ATMOS_RESTART_CHECK ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Data check : YES'
       ATMOS_sw_check = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Data check : NO'
       ATMOS_sw_check = .false.
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

    if ( ATMOS_QTRC_NMAX >= 1 ) then 

       A_VA = 5 + ATMOS_QTRC_NMAX
       allocate( A_NAME(A_VA) )
       allocate( A_DESC(A_VA) )
       allocate( A_UNIT(A_VA) )

       A_NAME( 1) = 'DENS'
       A_NAME( 2) = 'MOMX'
       A_NAME( 3) = 'MOMY'
       A_NAME( 4) = 'MOMZ'
       A_NAME( 5) = 'POTT'

       A_DESC( 1) = 'density'
       A_DESC( 2) = 'momentum (x)'
       A_DESC( 3) = 'momentum (y)'
       A_DESC( 4) = 'momentum (z)'
       A_DESC( 5) = 'potential temp.'

       A_UNIT( 1) = 'kg/m3'
       A_UNIT( 2) = 'kg/m2/s'
       A_UNIT( 3) = 'kg/m2/s'
       A_UNIT( 4) = 'kg/m2/s'
       A_UNIT( 5) = 'K'

       A_QA = ATMOS_QTRC_NMAX
       A_NAME(6:5+A_QA) = ATMOS_QTRC_VARNAME(1:ATMOS_QTRC_NMAX)
       A_DESC(6:5+A_QA) = ATMOS_QTRC_VARDESC(1:ATMOS_QTRC_NMAX)
       A_UNIT(6:5+A_QA) = ATMOS_QTRC_VARUNIT(1:ATMOS_QTRC_NMAX)

       if ( ATMOS_TYPE_PHY_MP == 'NDW6' ) then
          if ( ATMOS_QTRC_NMAX < 11 ) then
             write(*,*) 'xxx The number of tracers is not enough!', ATMOS_QTRC_NMAX
             call PRC_MPIstop
          endif

          A_NAME( 6) = 'QV'
          A_NAME( 7) = 'QC'
          A_NAME( 8) = 'QR'
          A_NAME( 9) = 'QI'
          A_NAME(10) = 'QS'
          A_NAME(11) = 'QG'
          A_NAME(12) = 'NC'
          A_NAME(13) = 'NR'
          A_NAME(14) = 'NI'
          A_NAME(15) = 'NS'
          A_NAME(16) = 'NG'

          A_DESC( 6) = 'Water Vapor mixing ratio'
          A_DESC( 7) = 'Cloud Water mixing ratio'
          A_DESC( 8) = 'Rain Water mixing ratio'
          A_DESC( 9) = 'Cloud Ice mixing ratio'
          A_DESC(10) = 'Snow mixing ratio'
          A_DESC(11) = 'Graupel mixing ratio'
          A_DESC(12) = 'Cloud Water Number Density'
          A_DESC(13) = 'Rain Water Number Density'
          A_DESC(14) = 'Cloud Ice Number Density'
          A_DESC(15) = 'Snow Number Density'
          A_DESC(16) = 'Graupel Number Density'

          A_UNIT( 6:11) = 'kg/kg'
          A_UNIT(12:16) = '1/m3'
       endif
    else
       write(*,*) 'xxx The number of tracers is not enough!', ATMOS_QTRC_NMAX
       A_QA = 0
    endif

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** [ATMOS] prognostic variables'
    if( IO_L ) write(IO_FID_LOG,*) &
    '***                 : VARNAME         , ', &
    'DESCRIPTION                                                     [UNIT            ]'
    do iv = 1, A_VA
       if( IO_L ) write(IO_FID_LOG,*) '*** NO.',iv,": ",A_NAME(iv),", ", A_DESC(iv),"[", A_UNIT(iv),"]"
    enddo

    allocate( atmos_var    (IA,JA,KA,A_VA) ); atmos_var(:,:,:,:)     = CONST_UNDEF8
    allocate( atmos_diagvar(IA,JA,KA,5)    ); atmos_diagvar(:,:,:,:) = CONST_UNDEF8

    I_QV = ATMOS_vars_getid( 'QV' )
    A_QWA = A_QWA + 1
    A_QWS = I_QV
    A_QWE = I_QV
    A_CPw(A_QWA) = CPvap
    A_CVw(A_QWA) = CVvap

    I_QC = ATMOS_vars_getid( 'QC' )
    if ( I_QC > 0 ) then
       A_QWA = A_QWA + 1
       A_QWE = I_QC
       A_QLS = I_QC
       A_QLE = I_QC
       A_CPw(A_QWA) = CL
       A_CVw(A_QWA) = CL
    endif
    I_QR = ATMOS_vars_getid( 'QR' )
    if ( I_QR > 0 ) then
       A_QWA = A_QWA + 1
       A_QWE = I_QR
       A_QLE = I_QR
       A_CPw(A_QWA) = CL
       A_CVw(A_QWA) = CL
    endif

    I_QI = ATMOS_vars_getid( 'QI' )
    if ( I_QI > 0 ) then
       A_QWA = A_QWA + 1
       A_QWE = I_QI
       A_QSS = I_QI
       A_QSE = I_QI
       A_CPw(A_QWA) = CI
       A_CVw(A_QWA) = CI
    endif
    I_QS = ATMOS_vars_getid( 'QS' )
    if ( I_QS > 0 ) then
       A_QWA = A_QWA + 1
       A_QWE = I_QS
       A_QSE = I_QS
       A_CPw(A_QWA) = CI
       A_CVw(A_QWA) = CI
    endif
    I_QG = ATMOS_vars_getid( 'QG' )
    if ( I_QG > 0 ) then
       A_QWA = A_QWA + 1
       A_QWE = I_QG
       A_QSE = I_QG
       A_CPw(A_QWA) = CI
       A_CVw(A_QWA) = CI
    endif

    I_NC = ATMOS_vars_getid( 'NC' )
    A_NWS = I_NC
    A_NWE = I_NC
    I_NR = ATMOS_vars_getid( 'NR' )
    if ( I_NR > 0 ) then
       A_NWE = I_NR
    endif
    I_NI = ATMOS_vars_getid( 'NI' )
    if ( I_NI > 0 ) then
       A_NWE = I_NI
    endif
    I_NS = ATMOS_vars_getid( 'NS' )
    if ( I_NS > 0 ) then
       A_NWE = I_NS
    endif
    I_NG = ATMOS_vars_getid( 'NG' )
    if ( I_NG > 0 ) then
       A_NWE = I_NG
    endif

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Hydrometeors Index Check'
    if( IO_L ) write(IO_FID_LOG,*) '*** Number of Hydrometeor tracers(QWA):', A_QWA
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I3,A,I3)') '*** Total  family QWS - QWE = ', A_QWS, ' - ', A_QWE
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I3,A,I3)') '*** Liquid family QLS - QLE = ', A_QLS, ' - ', A_QLE
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I3,A,I3)') '*** Soild  family QSS - QSE = ', A_QSS, ' - ', A_QSE
    if( IO_L ) write(IO_FID_LOG,*) '***  I_QV:', I_QV
    if( IO_L ) write(IO_FID_LOG,*) '***  I_QC:', I_QC
    if( IO_L ) write(IO_FID_LOG,*) '***  I_QR:', I_QR
    if( IO_L ) write(IO_FID_LOG,*) '***  I_QI:', I_QI
    if( IO_L ) write(IO_FID_LOG,*) '***  I_QS:', I_QS
    if( IO_L ) write(IO_FID_LOG,*) '***  I_QG:', I_QG
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I3,A,I3)') '*** Number family NWS -NWE = ', A_QWS, ' - ', A_QWE
    if( IO_L ) write(IO_FID_LOG,*) '***  I_NC:', I_NC
    if( IO_L ) write(IO_FID_LOG,*) '***  I_NR:', I_NR
    if( IO_L ) write(IO_FID_LOG,*) '***  I_NI:', I_NI
    if( IO_L ) write(IO_FID_LOG,*) '***  I_NS:', I_NS
    if( IO_L ) write(IO_FID_LOG,*) '***  I_NG:', I_NG

    return
  end subroutine ATMOS_vars_setup

  !-----------------------------------------------------------------------------
  !> search & get available tracer id by name
  !> @return fid
  !-----------------------------------------------------------------------------
  function ATMOS_vars_getid( name ) result(iv)
    implicit none

    character(len=*), intent(in) :: name !< tracer name
    integer                      :: iv   !< tracer ID
    !---------------------------------------------------------------------------

    do iv = 1, A_VA-5
       if( trim(A_NAME(iv+5)) == trim(name) ) return
    enddo

    iv = -1

  end function ATMOS_vars_getid

  !-----------------------------------------------------------------------------
  !> Read restart of atmospheric variables
  !-----------------------------------------------------------------------------
  subroutine ATMOS_vars_restart_read
    use mod_grid, only : &
       IA   => GRID_IA,   &
       JA   => GRID_JA,   &
       KA   => GRID_KA,   &
       IMAX => GRID_IMAX, &
       JMAX => GRID_JMAX, &
       KMAX => GRID_KMAX, &
       IS   => GRID_IS,   &
       IE   => GRID_IE,   &
       JS   => GRID_JS,   &
       JE   => GRID_JE,   &
       KS   => GRID_KS,   &
       KE   => GRID_KE
    use mod_comm, only: &
       COMM_vars,  &
       COMM_stats, &
       COMM_total
    use mod_fileio, only: &
       FIO_input
    implicit none

    real(8), allocatable :: restart_atmos(:,:,:) !> restart file (no HALO)

    character(len=IO_FILECHR) :: bname
    character(len=8)          :: lname

    integer :: iv, i, j
    !---------------------------------------------------------------------------

    allocate( restart_atmos(IMAX,JMAX,KMAX) )

    bname = ATMOS_RESTART_IN_BASENAME
    write(lname,'(A,I4.4)') 'ZDEF', KMAX

    do iv = 1, A_VA
       call FIO_input( restart_atmos(:,:,:), bname, A_NAME(iv), lname, 1, KMAX, 1 )

       atmos_var(IS:IE,JS:JE,KS:KE,iv) = restart_atmos(1:IMAX,1:JMAX,1:KMAX)
    enddo

    deallocate( restart_atmos )

    ! fill IHALO & JHALO
    call COMM_vars( atmos_var(:,:,:,:) )

    ! fill KHALO
    do iv = 1, A_VA
    do j  = 1, JA
    do i  = 1, IA
       atmos_var(i,j,   1:KS-1,iv) = atmos_var(i,j,KS,iv)
       atmos_var(i,j,KE+1:KA,  iv) = atmos_var(i,j,KE,iv)
    enddo
    enddo
    enddo

    call COMM_stats( atmos_var(:,:,:,:), A_NAME(:) )

    ! check total mass
    call COMM_total( atmos_var(:,:,:,1:1), A_NAME(1:1) )

    ! atmos_var -> atmos_diagvar
    call ATMOS_DMP2PVT

    ! fill IHALO & JHALO
    call COMM_vars( atmos_diagvar(:,:,:,:) )

    return
  end subroutine ATMOS_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Write restart of atmospheric variables
  !-----------------------------------------------------------------------------
  subroutine ATMOS_vars_restart_write
    use mod_stdio, only: &
       IO_FID_LOG,  &
       IO_L
    use mod_time, only: &
       NOWSEC => TIME_NOWSEC
    use mod_grid, only : &
       IMAX => GRID_IMAX, &
       JMAX => GRID_JMAX, &
       KMAX => GRID_KMAX, &
       IS   => GRID_IS,   &
       IE   => GRID_IE,   &
       JS   => GRID_JS,   &
       JE   => GRID_JE,   &
       KS   => GRID_KS,   &
       KE   => GRID_KE
    use mod_comm, only: &
       COMM_vars,  &
       COMM_stats, &
       COMM_total
    use mod_fileio, only: &
       FIO_output
    implicit none

    real(8), allocatable :: restart_atmos(:,:,:) !> restart file (no HALO)

    character(len=IO_FILECHR) :: bname
    character(len=FIO_HMID)   :: desc
    character(len=8)          :: lname

    integer :: iv
    !---------------------------------------------------------------------------

    allocate( restart_atmos(IMAX,JMAX,KMAX) )

    call COMM_stats( atmos_var(:,:,:,:), A_NAME(:) )

    ! check total mass
    call COMM_total( atmos_var(:,:,:,1:1), A_NAME(1:1) )

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Output restart file (atmos) ***'

    write(bname,'(A,A,F15.3)') trim(ATMOS_RESTART_OUT_BASENAME), '_', NOWSEC
    desc  = 'SCALE3 PROGNOSTIC VARS.'
    write(lname,'(A,I4.4)') 'ZDEF', KMAX

    do iv = 1, A_VA
       restart_atmos(1:IMAX,1:JMAX,1:KMAX) = atmos_var(IS:IE,JS:JE,KS:KE,iv)

       call FIO_output( restart_atmos(:,:,:), bname, desc, '',       &
                        A_NAME(iv), A_DESC(iv), '', A_UNIT(iv),      &
                        FIO_REAL8, lname, 1, KMAX, 1, NOWSEC, NOWSEC )
    enddo

    deallocate( restart_atmos )

    return
  end subroutine ATMOS_vars_restart_write

  !-----------------------------------------------------------------------------
  !> Check and compare between last data and sample data
  !-----------------------------------------------------------------------------
  subroutine ATMOS_vars_restart_check
    use mod_stdio, only: &
       IO_FID_LOG,  &
       IO_L
    use mod_process, only: &
       PRC_myrank
    use mod_grid, only : &
       IA   => GRID_IA,   &
       JA   => GRID_JA,   &
       KA   => GRID_KA,   &
       IMAX => GRID_IMAX, &
       JMAX => GRID_JMAX, &
       KMAX => GRID_KMAX, &
       IS   => GRID_IS,   &
       IE   => GRID_IE,   &
       JS   => GRID_JS,   &
       JE   => GRID_JE,   &
       KS   => GRID_KS,   &
       KE   => GRID_KE
    use mod_comm, only: &
       COMM_total
    use mod_fileio, only: &
       FIO_input
    implicit none

    real(8), allocatable :: restart_atmos(:,:,:) !> restart file (no HALO)

    real(8), allocatable :: atmos_var_check(:,:,:,:)

    character(len=IO_FILECHR) :: bname
    character(len=FIO_HMID)   :: desc
    character(len=8)          :: lname

    logical :: datacheck
    integer :: i, j, k, iv
    !---------------------------------------------------------------------------

    allocate( restart_atmos(IMAX,JMAX,KMAX) )
    allocate( atmos_var_check(IA,JA,KA,A_VA) )

    bname = ATMOS_RESTART_CHECK_BASENAME
    write(lname,'(A,I4.4)') 'ZDEF', KMAX

    do iv = 1, A_VA
       call FIO_input( restart_atmos(:,:,:), bname, A_NAME(iv), lname, 1, KMAX, 1 )

       atmos_var_check(IS:IE,JS:JE,KS:KE,iv) = restart_atmos(1:IMAX,1:JMAX,1:KMAX)
    enddo

    deallocate( restart_atmos )

    ! check total mass
    call COMM_total( atmos_var(:,:,:,1:1), A_NAME(1:1) )

    write(*,*) 'Compare last Data with ', trim(ATMOS_RESTART_CHECK_BASENAME), 'on PE=', PRC_myrank
    write(*,*) '*** criterion = ', ATMOS_RESTART_CHECK_CRITERION
    datacheck = .true.

    do iv = 1, A_VA
    do k = KS, KE
    do j = JS, JE
    do i = IS, IE
       if ( abs( atmos_var(i,j,k,iv)-atmos_var_check(i,j,k,iv) ) > ATMOS_RESTART_CHECK_CRITERION ) then
          write(*,*) 'xxx there is the difference : ', atmos_var(i,j,k,iv)-atmos_var_check(i,j,k,iv)
          write(*,*) 'xxx at (PE-id,i,j,k,varname) : ', PRC_myrank, i, j, k, A_NAME(iv)
          datacheck = .false.
       endif
    enddo
    enddo
    enddo
    enddo
    if (datacheck) then
       if( IO_L ) write(IO_FID_LOG,*) 'Data Check Clear.'
       write(*,*) 'Data Check Clear.'
    endif
    call COMM_total( atmos_var_check(:,:,:,1:1), A_NAME(1:1) )

    return
  end subroutine ATMOS_vars_restart_check

  !-----------------------------------------------------------------------------
  !> Put and Communicate prognostic variables
  !-----------------------------------------------------------------------------
  subroutine ATMOS_vars_put( &
       dens, &
       momx, &
       momy, &
       momz, &
       pott, &
       qtrc  )
    use mod_grid, only: &
       IA   => GRID_IA,   &
       JA   => GRID_JA,   &
       KA   => GRID_KA
    use mod_comm, only: &
       COMM_vars
    implicit none

    real(8), intent(in) :: dens(IA,JA,KA)
    real(8), intent(in) :: momx(IA,JA,KA)
    real(8), intent(in) :: momy(IA,JA,KA)
    real(8), intent(in) :: momz(IA,JA,KA)
    real(8), intent(in) :: pott(IA,JA,KA)

    real(8), intent(in) :: qtrc(IA,JA,KA,A_QA)

    integer :: iq
    !---------------------------------------------------------------------------

    atmos_var(:,:,:,1) = dens(:,:,:)
    atmos_var(:,:,:,2) = momx(:,:,:)
    atmos_var(:,:,:,3) = momy(:,:,:)
    atmos_var(:,:,:,4) = momz(:,:,:)
    atmos_var(:,:,:,5) = pott(:,:,:)

    if ( A_QA > 0 ) then
       do iq = 1, A_QA
          atmos_var(:,:,:,5+iq) = qtrc(:,:,:,iq)
       enddo
    endif

    ! fill IHALO & JHALO
    call COMM_vars( atmos_var(:,:,:,:) )

    ! atmos_var -> atmos_diagvar
    call ATMOS_DMP2PVT

    ! fill IHALO & JHALO
    call COMM_vars( atmos_diagvar(:,:,:,:) )

    return
  end subroutine ATMOS_vars_put

  !-----------------------------------------------------------------------------
  !> Get prognostic variables
  !-----------------------------------------------------------------------------
  subroutine ATMOS_vars_get( &
       dens, &
       momx, &
       momy, &
       momz, &
       pott, &
       qtrc  )
    use mod_grid, only: &
       IA => GRID_IA, &
       JA => GRID_JA, &
       KA => GRID_KA
    implicit none

    real(8), intent(out) :: dens(IA,JA,KA)
    real(8), intent(out) :: momx(IA,JA,KA)
    real(8), intent(out) :: momy(IA,JA,KA)
    real(8), intent(out) :: momz(IA,JA,KA)
    real(8), intent(out) :: pott(IA,JA,KA)

    real(8), intent(out) :: qtrc(IA,JA,KA,A_QA)

    integer :: iq
    !---------------------------------------------------------------------------

    dens(:,:,:) = atmos_var(:,:,:,1)
    momx(:,:,:) = atmos_var(:,:,:,2)
    momy(:,:,:) = atmos_var(:,:,:,3)
    momz(:,:,:) = atmos_var(:,:,:,4)
    pott(:,:,:) = atmos_var(:,:,:,5)

    if ( A_QA > 0 ) then
       do iq = 1, A_QA
          qtrc(:,:,:,iq) = atmos_var(:,:,:,5+iq)
       enddo
    endif

    return
  end subroutine ATMOS_vars_get

  !-----------------------------------------------------------------------------
  !> Get prognostic variables
  !-----------------------------------------------------------------------------
  subroutine ATMOS_vars_getall( &
       dens, &
       momx, &
       momy, &
       momz, &
       pott, &
       qtrc, &
       pres, &
       velx, &
       vely, &
       velz, &
       temp  )
    use mod_grid, only: &
       IA => GRID_IA, &
       JA => GRID_JA, &
       KA => GRID_KA
    implicit none

    real(8), intent(out) :: dens(IA,JA,KA)
    real(8), intent(out) :: momx(IA,JA,KA)
    real(8), intent(out) :: momy(IA,JA,KA)
    real(8), intent(out) :: momz(IA,JA,KA)
    real(8), intent(out) :: pott(IA,JA,KA)

    real(8), intent(out) :: qtrc(IA,JA,KA,A_QA)

    real(8), intent(out) :: pres(IA,JA,KA)
    real(8), intent(out) :: velx(IA,JA,KA)
    real(8), intent(out) :: vely(IA,JA,KA)
    real(8), intent(out) :: velz(IA,JA,KA)
    real(8), intent(out) :: temp(IA,JA,KA)

    integer :: iq
    !---------------------------------------------------------------------------

    dens(:,:,:) = atmos_var(:,:,:,1)
    momx(:,:,:) = atmos_var(:,:,:,2)
    momy(:,:,:) = atmos_var(:,:,:,3)
    momz(:,:,:) = atmos_var(:,:,:,4)
    pott(:,:,:) = atmos_var(:,:,:,5)

    if ( A_QA > 0 ) then
       do iq = 1, A_QA
          qtrc(:,:,:,iq) = atmos_var(:,:,:,5+iq)
       enddo
    endif

    pres(:,:,:) = atmos_diagvar(:,:,:,1)
    velx(:,:,:) = atmos_diagvar(:,:,:,2)
    vely(:,:,:) = atmos_diagvar(:,:,:,3)
    velz(:,:,:) = atmos_diagvar(:,:,:,4)
    temp(:,:,:) = atmos_diagvar(:,:,:,5)

    return
  end subroutine ATMOS_vars_getall

  !-----------------------------------------------------------------------------
  !> Get prognostic variables
  !-----------------------------------------------------------------------------
  subroutine ATMOS_DMP2PVT
    use mod_stdio, only: &
       IO_FID_LOG,  &
       IO_L
    use mod_const, only : &
       Rair   => CONST_Rair,   &
       CPair  => CONST_CPair,  &
       RovCP  => CONST_RovCP,  &
       CPovCV => CONST_CPovCV, &
       LH0    => CONST_LH0,    &
       Pstd   => CONST_Pstd
    use mod_grid, only: &
       IA   => GRID_IA,   &
       JA   => GRID_JA,   &
       KA   => GRID_KA,   &
       IS   => GRID_IS,   &
       IE   => GRID_IE,   &
       JS   => GRID_JS,   &
       JE   => GRID_JE,   &
       KS   => GRID_KS,   &
       KE   => GRID_KE,   &
       WS   => GRID_WS,   &
       WE   => GRID_WE
    implicit none

    real(8) :: pres(IA,JA,KA)
    real(8) :: velx(IA,JA,KA)
    real(8) :: vely(IA,JA,KA)
    real(8) :: velz(IA,JA,KA)
    real(8) :: temp(IA,JA,KA)

    integer :: i, j, k
    !---------------------------------------------------------------------------

    ! momentum -> velocity
    do k = KS, KE
    do j = JS, JE
    do i = IS, IE
       velx(i,j,k) = 2.D0 * atmos_var(i,j,k,2) / ( atmos_var(i+1,j,  k,1)+atmos_var(i,j,k,1) )
       vely(i,j,k) = 2.D0 * atmos_var(i,j,k,3) / ( atmos_var(i,  j+1,k,1)+atmos_var(i,j,k,1) )
    enddo
    enddo
    enddo

    do k = WS+1, WE-1
    do j = JS,   JE
    do i = IS,   IE
       velz(i,j,k) = 2.D0 * atmos_var(i,j,k,4) / ( atmos_var(i,j,k+1,1)+atmos_var(i,j,k,1) )
    enddo
    enddo
    enddo
    velz(:,:,WS) = 0.D0 ! bottom boundary
    velz(:,:,WE) = 0.D0 ! top    boundary

    ! diagnose pressure, temperature
    do k = KS, KE
    do j = JS, JE
    do i = IS, IE
       pres(i,j,k) = Pstd * ( atmos_var(i,j,k,1) * atmos_var(i,j,k,5) * Rair / Pstd )**CPovCV ! first guess

       temp(i,j,k) = pres(i,j,k) / ( atmos_var(i,j,k,1) * Rair )
    enddo
    enddo
    enddo

    atmos_diagvar(IS:IE,JS:JE,KS:KE,1) = pres(IS:IE,JS:JE,KS:KE)
    atmos_diagvar(IS:IE,JS:JE,KS:KE,2) = velx(IS:IE,JS:JE,KS:KE)
    atmos_diagvar(IS:IE,JS:JE,KS:KE,3) = vely(IS:IE,JS:JE,KS:KE)
    atmos_diagvar(IS:IE,JS:JE,WS:WE,4) = velz(IS:IE,JS:JE,WS:WE)
    atmos_diagvar(IS:IE,JS:JE,KS:KE,5) = temp(IS:IE,JS:JE,KS:KE)

    return
  end subroutine ATMOS_DMP2PVT

end module mod_atmos_vars
