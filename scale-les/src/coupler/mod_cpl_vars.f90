!-------------------------------------------------------------------------------
!> module COUPLER Variables
!!
!! @par Description
!!          Container for coupler variables
!!
!! @author Team SCALE
!! @li      2013-08-31 (T.Yamaura)  [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_cpl_vars
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
  public :: CPL_vars_setup
  public :: CPL_vars_fillhalo
  public :: CPL_vars_restart_read
  public :: CPL_vars_restart_write
  public :: CPL_vars_history
  public :: CPL_vars_total

  public :: CPL_vars_merge
  public :: CPL_putAtm
  public :: CPL_putLnd
  public :: CPL_putUrb
  public :: CPL_putOcn
  public :: CPL_getAtm
  public :: CPL_getAtm_RD
  public :: CPL_getLnd
  public :: CPL_getUrb
  public :: CPL_getOcn

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical,                public, save :: CPL_do          = .true.  ! main switch for the model

  character(len=H_SHORT), public, save :: CPL_TYPE_AtmLnd = 'OFF'   !< atmos-land coupler type
  character(len=H_SHORT), public, save :: CPL_TYPE_AtmUrb = 'OFF'   !< atmos-urban coupler type
  character(len=H_SHORT), public, save :: CPL_TYPE_AtmOcn = 'OFF'   !< atmos-ocean coupler type
  logical,                public, save :: CPL_LST_UPDATE  = .false. !< is Land Surface Temperature updated?
  logical,                public, save :: CPL_UST_UPDATE  = .false. !< is Urban Surface Temperature updated?
  logical,                public, save :: CPL_SST_UPDATE  = .false. !< is Sea  Surface Temperature updated?

  logical,                public, save :: CPL_sw_ALL      = .false. !< do coupler calculation?
  logical,                public, save :: CPL_sw_AtmOcn             !< do atmos-ocean coupler calculation?
  logical,                public, save :: CPL_sw_AtmLnd             !< do atmos-land  coupler calculation?
  logical,                public, save :: CPL_sw_AtmUrb             !< do atmos-urban coupler calculation?
  logical,                public, save :: CPL_sw_restart            !< output coupler restart?

  ! prognostic variables
  real(RP), public, save, allocatable :: LST  (:,:)   ! Land Surface Temperature [K]
  real(RP), public, save, allocatable :: UST  (:,:)   ! Urban Surface Temperature [K]
  real(RP), public, save, allocatable :: SST  (:,:)   ! Sea Surface Temperature [K]
  real(RP), public, save, allocatable :: ALBG (:,:,:) ! Land Surface albedo [0-1]
  real(RP), public, save, allocatable :: ALBW (:,:,:) ! Sea Surface albedo [0-1]
  real(RP), public, save, allocatable :: Z0W  (:,:)   ! Sea Surface roughness length [m]
  real(RP), public, save, allocatable :: SkinT(:,:)   ! Ground Skin Temperature [K]
  real(RP), public, save, allocatable :: SkinW(:,:)   ! Ground Skin Water       [kg/m2]
  real(RP), public, save, allocatable :: SnowQ(:,:)   ! Ground Snow amount      [kg/m2]
  real(RP), public, save, allocatable :: SnowT(:,:)   ! Ground Snow Temperature [K]

  ! surface fluxes to atmospheric model
  real(RP), public, save, allocatable :: Atm_XMFLX (:,:) ! x-momentum flux [kg/m2/s]
  real(RP), public, save, allocatable :: Atm_YMFLX (:,:) ! y-momentum flux [kg/m2/s]
  real(RP), public, save, allocatable :: Atm_ZMFLX (:,:) ! z-momentum flux [kg/m2/s]
  real(RP), public, save, allocatable :: Atm_SHFLX (:,:) ! sensible heat flux (upward positive) [W/m2]
  real(RP), public, save, allocatable :: Atm_LHFLX (:,:) ! latent heat flux (upward positive) [W/m2]
  real(RP), public, save, allocatable :: Atm_QVFLX (:,:) ! moisture flux for atmosphere [kg/m2/s]

  ! surface fluxes to land model
  real(RP), public, save, allocatable :: Lnd_GHFLX  (:,:) ! ground heat flux (upward positive) [W/m2]
  real(RP), public, save, allocatable :: Lnd_PRECFLX(:,:) ! precipitation flux [kg/m2/s]
  real(RP), public, save, allocatable :: Lnd_QVFLX  (:,:) ! moisture flux for land [kg/m2/s]

  ! surface fluxes to urban model
  real(RP), public, save, allocatable :: Urb_GHFLX  (:,:) ! ground heat flux (upward positive) [W/m2]
  real(RP), public, save, allocatable :: Urb_PRECFLX(:,:) ! precipitation flux [kg/m2/s]
  real(RP), public, save, allocatable :: Urb_QVFLX  (:,:) ! moisture flux for land [kg/m2/s]

  ! surface fluxes to ocean model
  real(RP), public, save, allocatable :: Ocn_WHFLX  (:,:) ! water heat flux (upward positive) [W/m2]
  real(RP), public, save, allocatable :: Ocn_PRECFLX(:,:) ! precipitation flux [kg/m2/s]
  real(RP), public, save, allocatable :: Ocn_QVFLX  (:,:) ! moisture flux for ocean [kg/m2/s]

  ! surface fluxes from atmosphere-ocean coupler
  real(RP), public, save, allocatable :: CPL_AtmOcn_XMFLX (:,:) ! x-momentum flux [kg/m2/s]
  real(RP), public, save, allocatable :: CPL_AtmOcn_YMFLX (:,:) ! y-momentum flux [kg/m2/s]
  real(RP), public, save, allocatable :: CPL_AtmOcn_ZMFLX (:,:) ! z-momentum flux [kg/m2/s]
  real(RP), public, save, allocatable :: CPL_AtmOcn_SHFLX (:,:) ! sensible heat flux (upward positive) [W/m2]
  real(RP), public, save, allocatable :: CPL_AtmOcn_LHFLX (:,:) ! latent   heat flux (upward positive) [W/m2]
  real(RP), public, save, allocatable :: CPL_AtmOcn_QVFLX (:,:) ! moisture flux for atmosphere [kg/m2/s]

  ! surface fluxes from atmosphere-land coupler
  real(RP), public, save, allocatable :: CPL_AtmLnd_XMFLX (:,:) ! x-momentum flux [kg/m2/s]
  real(RP), public, save, allocatable :: CPL_AtmLnd_YMFLX (:,:) ! y-momentum flux [kg/m2/s]
  real(RP), public, save, allocatable :: CPL_AtmLnd_ZMFLX (:,:) ! z-momentum flux [kg/m2/s]
  real(RP), public, save, allocatable :: CPL_AtmLnd_SHFLX (:,:) ! sensible heat flux (upward positive) [W/m2]
  real(RP), public, save, allocatable :: CPL_AtmLnd_LHFLX (:,:) ! latent   heat flux (upward positive) [W/m2]
  real(RP), public, save, allocatable :: CPL_AtmLnd_QVFLX (:,:) ! moisture flux for atmosphere [kg/m2/s]

  ! surface fluxes from atmosphere-urban coupler
  real(RP), public, save, allocatable :: CPL_AtmUrb_XMFLX (:,:) ! x-momentum flux [kg/m2/s]
  real(RP), public, save, allocatable :: CPL_AtmUrb_YMFLX (:,:) ! y-momentum flux [kg/m2/s]
  real(RP), public, save, allocatable :: CPL_AtmUrb_ZMFLX (:,:) ! z-momentum flux [kg/m2/s]
  real(RP), public, save, allocatable :: CPL_AtmUrb_SHFLX (:,:) ! sensible heat flux (upward positive) [W/m2]
  real(RP), public, save, allocatable :: CPL_AtmUrb_LHFLX (:,:) ! latent   heat flux (upward positive) [W/m2]
  real(RP), public, save, allocatable :: CPL_AtmUrb_QVFLX (:,:) ! moisture flux for atmosphere [kg/m2/s]

  ! Atmospheric values
  real(RP), public, save, allocatable :: CPL_DENS(:,:) ! air density [kg/m3]
  real(RP), public, save, allocatable :: CPL_MOMX(:,:) ! momentum x [kg/m2/s]
  real(RP), public, save, allocatable :: CPL_MOMY(:,:) ! momentum y [kg/m2/s]
  real(RP), public, save, allocatable :: CPL_MOMZ(:,:) ! momentum z [kg/m2/s]
  real(RP), public, save, allocatable :: CPL_TMPA(:,:) ! temperature at the lowest atmospheric layer [K]
  real(RP), public, save, allocatable :: CPL_PRSA(:,:) ! pressure at the lowest atmospheric layer [Pa]
  real(RP), public, save, allocatable :: CPL_QV  (:,:) ! ratio of mass of tracer to total mass [kg/kg]
  real(RP), public, save, allocatable :: CPL_PRSS(:,:) ! pressure at the surface [Pa]
  real(RP), public, save, allocatable :: CPL_PREC(:,:) ! surface precipitation rate [kg/m2/s]
  real(RP), public, save, allocatable :: CPL_SWD (:,:) ! downward short-wave radiation flux (upward positive) [W/m2]
  real(RP), public, save, allocatable :: CPL_LWD (:,:) ! downward long-wave radiation flux (upward positive) [W/m2]

  ! Ocean values
  real(RP), public, save, allocatable :: CPL_TW  (:,:) ! water temperature [K]

  ! Land values
  real(RP), public, save, allocatable :: CPL_TG  (:,:) ! soil temperature [K]
  real(RP), public, save, allocatable :: CPL_QVEF(:,:) ! efficiency of evaporation [0-1]
  real(RP), public, save, allocatable :: CPL_TCS (:,:) ! thermal conductivity for soil [W/m/K]
  real(RP), public, save, allocatable :: CPL_DZG (:,:) ! soil depth [m]
  real(RP), public, save, allocatable :: CPL_Z0M (:,:) ! roughness length for momemtum [m]
  real(RP), public, save, allocatable :: CPL_Z0H (:,:) ! roughness length for heat [m]
  real(RP), public, save, allocatable :: CPL_Z0E (:,:) ! roughness length for vapor [m]

  ! counter
  real(RP), public, save :: CNT_Atm_Lnd ! counter for atmos flux by land
  real(RP), public, save :: CNT_Atm_Urb ! counter for atmos flux by urban
  real(RP), public, save :: CNT_Atm_Ocn ! counter for atmos flux by ocean
  real(RP), public, save :: CNT_Lnd     ! counter for land flux
  real(RP), public, save :: CNT_Urb     ! counter for urban flux
  real(RP), public, save :: CNT_Ocn     ! counter for ocean flux

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical,                private, save      :: CPL_RESTART_OUTPUT       = .false.       !< output restart file?
  character(len=H_LONG),  private, save      :: CPL_RESTART_IN_BASENAME  = ''            !< basename of the input file
  character(len=H_LONG),  private, save      :: CPL_RESTART_OUT_BASENAME = ''            !< basename of the output file
  character(len=H_MID),   private, save      :: CPL_RESTART_OUT_TITLE    = 'CPL restart' !< title    of the output file
  character(len=H_MID),   private, save      :: CPL_RESTART_OUT_DTYPE    = 'DEFAULT'     !< REAL4 or REAL8

  logical,                private, save      :: CPL_VARS_CHECKRANGE      = .false.

  integer,                private, parameter :: VMAX          = 12
  integer,                private, parameter :: I_LST         =  1
  integer,                private, parameter :: I_UST         =  2
  integer,                private, parameter :: I_SST         =  3
  integer,                private, parameter :: I_ALBW_SW     =  4
  integer,                private, parameter :: I_ALBW_LW     =  5
  integer,                private, parameter :: I_ALBG_SW     =  6
  integer,                private, parameter :: I_ALBG_LW     =  7
  integer,                private, parameter :: I_Z0W         =  8
  integer,                private, parameter :: I_SkinT       =  9
  integer,                private, parameter :: I_SkinW       = 10
  integer,                private, parameter :: I_SnowQ       = 11
  integer,                private, parameter :: I_SnowT       = 12

  character(len=H_SHORT), private, save      :: VAR_NAME(VMAX) !< name  of the coupler variables
  character(len=H_MID),   private, save      :: VAR_DESC(VMAX) !< desc. of the coupler variables
  character(len=H_SHORT), private, save      :: VAR_UNIT(VMAX) !< unit  of the coupler variables

  data VAR_NAME / 'LST',         &
                  'UST',         &
                  'SST',         &
                  'ALBW_SW',     &
                  'ALBW_LW',     &
                  'ALBG_SW',     &
                  'ALBG_LW',     &
                  'Z0W',         &
                  'SkinT',       &
                  'SkinW',       &
                  'SnowQ',       &
                  'SnowT'        /

  data VAR_DESC / 'land surface temp.',               &
                  'urban surface temp.',              &
                  'sea surface temp.',                &
                  'sea surface albedo for SW',        &
                  'sea surface albedo for LW',        &
                  'land surface albedo for SW',       &
                  'land surface albedo for LW',       &
                  'sea surface roughness length',     &
                  'ground skin temp.',                &
                  'ground skin water',                &
                  'ground snow amount',               &
                  'ground snow temp.'                 /

  data VAR_UNIT / 'K',       &
                  'K',       &
                  'K',       &
                  '0-1',     &
                  '0-1',     &
                  '0-1',     &
                  '0-1',     &
                  'm',       &
                  'K',       &
                  'kg/m2',   &
                  'kg/m2',   &
                  'K'        /

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine CPL_vars_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    NAMELIST / PARAM_CPL / &
       CPL_do,          &
       CPL_TYPE_AtmLnd, &
       CPL_TYPE_AtmUrb, &
       CPL_TYPE_AtmOcn, &
       CPL_LST_UPDATE,  &
       CPL_UST_UPDATE,  &
       CPL_SST_UPDATE

    NAMELIST / PARAM_CPL_VARS /  &
       CPL_RESTART_IN_BASENAME,  &
       CPL_RESTART_OUTPUT,       &
       CPL_RESTART_OUT_BASENAME, &
       CPL_RESTART_OUT_TITLE,    &
       CPL_RESTART_OUT_DTYPE,    &
       CPL_VARS_CHECKRANGE

    integer :: ip, ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[CPL VARS]/Categ[CPL]'
    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_CPL,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_CPL. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_CPL)

    if( IO_L ) write(IO_FID_LOG,*) '*** [CPL] selected components'

    ! Atoms-Land Switch
    if ( CPL_TYPE_AtmLnd /= 'OFF' .AND. CPL_TYPE_AtmLnd /= 'NONE' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Atmos-Land Coupler : ON'
       CPL_sw_AtmLnd = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Atmos-Land Coupler : OFF'
       CPL_sw_AtmLnd = .false.
    endif

    ! Atoms-Urban Switch
    if ( CPL_TYPE_AtmUrb /= 'OFF' .AND. CPL_TYPE_AtmUrb /= 'NONE' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Atmos-Urban Coupler : ON'
       CPL_sw_AtmUrb = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Atmos-Urban Coupler : OFF'
       CPL_sw_AtmUrb = .false.
    endif

    ! Atoms-Ocean Switch
    if ( CPL_TYPE_AtmOcn /= 'OFF' .AND. CPL_TYPE_AtmOcn /= 'NONE' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Atmos-Ocean Coupler : ON'
       CPL_sw_AtmOcn = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Atmos-Ocean Coupler : OFF'
       CPL_sw_AtmOcn = .false.
    endif

    if ( CPL_sw_AtmLnd .OR. CPL_sw_AtmOcn ) then
       CPL_sw_ALL = .true.
    endif

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_CPL_VARS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_CPL_VARS. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_CPL_VARS)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** [CPL] prognostic variables'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A8,A,A32,3(A))') &
               '***       |',' VARNAME','|', 'DESCRIPTION                     ','[', 'UNIT            ',']'
    do ip = 1, VMAX
       if( IO_L ) write(IO_FID_LOG,'(1x,A,i3,A,A8,A,A32,3(A))') &
                  '*** NO.',ip,'|',trim(VAR_NAME(ip)),'|', VAR_DESC(ip),'[', VAR_UNIT(ip),']'
    enddo

    ! restart switch
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) 'Output...'
    if ( CPL_RESTART_OUTPUT ) then
       if( IO_L ) write(IO_FID_LOG,*) '  Coupler restart output : YES'
       CPL_sw_restart = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '  Coupler restart output : NO'
       CPL_sw_restart = .false.
    endif
    if( IO_L ) write(IO_FID_LOG,*)

    CNT_Atm_Lnd = 0.0_RP
    CNT_Atm_Urb = 0.0_RP
    CNT_Atm_Ocn = 0.0_RP
    CNT_Lnd     = 0.0_RP
    CNT_Urb     = 0.0_RP
    CNT_Ocn     = 0.0_RP

    allocate( LST  (IA,JA)   )
    allocate( UST  (IA,JA)   )
    allocate( SST  (IA,JA)   )
    allocate( ALBW (IA,JA,2) )
    allocate( ALBG (IA,JA,2) )
    allocate( Z0W  (IA,JA)   )
    allocate( SkinT(IA,JA)   )
    allocate( SkinW(IA,JA)   )
    allocate( SnowQ(IA,JA)   )
    allocate( SnowT(IA,JA)   )

    ! work variables
    allocate( Atm_XMFLX(IA,JA) )
    allocate( Atm_YMFLX(IA,JA) )
    allocate( Atm_ZMFLX(IA,JA) )
    allocate( Atm_SHFLX(IA,JA) )
    allocate( Atm_LHFLX(IA,JA) )
    allocate( Atm_QVFLX(IA,JA) )

    allocate( Lnd_GHFLX  (IA,JA) )
    allocate( Lnd_PRECFLX(IA,JA) )
    allocate( Lnd_QVFLX  (IA,JA) )

    allocate( Urb_GHFLX  (IA,JA) )
    allocate( Urb_PRECFLX(IA,JA) )
    allocate( Urb_QVFLX  (IA,JA) )

    allocate( Ocn_WHFLX  (IA,JA) )
    allocate( Ocn_PRECFLX(IA,JA) )
    allocate( Ocn_QVFLX  (IA,JA) )

    allocate( CPL_DENS(IA,JA) )
    allocate( CPL_MOMX(IA,JA) )
    allocate( CPL_MOMY(IA,JA) )
    allocate( CPL_MOMZ(IA,JA) )
    allocate( CPL_TMPA(IA,JA) )
    allocate( CPL_PRSA(IA,JA) )
    allocate( CPL_QV  (IA,JA) )
    allocate( CPL_PRSS(IA,JA) )
    allocate( CPL_PREC(IA,JA) )
    allocate( CPL_SWD (IA,JA) )
    allocate( CPL_LWD (IA,JA) )

    allocate( CPL_TG  (IA,JA) )
    allocate( CPL_QVEF(IA,JA) )
    allocate( CPL_TCS (IA,JA) )
    allocate( CPL_DZG (IA,JA) )
    allocate( CPL_Z0M (IA,JA) )
    allocate( CPL_Z0H (IA,JA) )
    allocate( CPL_Z0E (IA,JA) )

    allocate( CPL_TW  (IA,JA) )

    allocate( CPL_AtmLnd_XMFLX(IA,JA) )
    allocate( CPL_AtmLnd_YMFLX(IA,JA) )
    allocate( CPL_AtmLnd_ZMFLX(IA,JA) )
    allocate( CPL_AtmLnd_SHFLX(IA,JA) )
    allocate( CPL_AtmLnd_LHFLX(IA,JA) )
    allocate( CPL_AtmLnd_QVFLX(IA,JA) )

    allocate( CPL_AtmUrb_XMFLX(IA,JA) )
    allocate( CPL_AtmUrb_YMFLX(IA,JA) )
    allocate( CPL_AtmUrb_ZMFLX(IA,JA) )
    allocate( CPL_AtmUrb_SHFLX(IA,JA) )
    allocate( CPL_AtmUrb_LHFLX(IA,JA) )
    allocate( CPL_AtmUrb_QVFLX(IA,JA) )

    allocate( CPL_AtmOcn_XMFLX(IA,JA) )
    allocate( CPL_AtmOcn_YMFLX(IA,JA) )
    allocate( CPL_AtmOcn_ZMFLX(IA,JA) )
    allocate( CPL_AtmOcn_SHFLX(IA,JA) )
    allocate( CPL_AtmOcn_LHFLX(IA,JA) )
    allocate( CPL_AtmOcn_QVFLX(IA,JA) )

    return
  end subroutine CPL_vars_setup

  !-----------------------------------------------------------------------------
  !> HALO Communication
  subroutine CPL_vars_fillhalo
    use scale_const, only: &
       I_SW => CONST_I_SW, &
       I_LW => CONST_I_LW
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none
    !---------------------------------------------------------------------------

    call COMM_vars8( LST  (:,:),       1 )
    call COMM_vars8( UST  (:,:),       2 )
    call COMM_vars8( SST  (:,:),       3 )
    call COMM_vars8( ALBW (:,:,I_SW),  4 )
    call COMM_vars8( ALBW (:,:,I_LW),  5 )
    call COMM_vars8( ALBG (:,:,I_SW),  6 )
    call COMM_vars8( ALBG (:,:,I_LW),  7 )
    call COMM_vars8( Z0W  (:,:),       8 )
    call COMM_vars8( SkinT(:,:),       9 )
    call COMM_vars8( SkinW(:,:),      10 )
    call COMM_vars8( SnowQ(:,:),      11 )
    call COMM_vars8( SnowT(:,:),      12 )

    call COMM_wait ( LST  (:,:),       1 )
    call COMM_wait ( UST  (:,:),       2 )
    call COMM_wait ( SST  (:,:),       3 )
    call COMM_wait ( ALBW (:,:,I_SW),  4 )
    call COMM_wait ( ALBW (:,:,I_LW),  5 )
    call COMM_wait ( ALBG (:,:,I_SW),  6 )
    call COMM_wait ( ALBG (:,:,I_LW),  7 )
    call COMM_wait ( Z0W  (:,:),       8 )
    call COMM_wait ( SkinT(:,:),       9 )
    call COMM_wait ( SkinW(:,:),      10 )
    call COMM_wait ( SnowQ(:,:),      11 )
    call COMM_wait ( SnowT(:,:),      12 )

    return
  end subroutine CPL_vars_fillhalo

  !-----------------------------------------------------------------------------
  !> Read coupler restart
  subroutine CPL_vars_restart_read
    use scale_const, only: &
       UNDEF => CONST_UNDEF, &
       I_SW  => CONST_I_SW,  &
       I_LW  => CONST_I_LW
    use scale_fileio, only: &
       FILEIO_read
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input restart file (coupler) ***'

    if ( CPL_RESTART_IN_BASENAME /= '' ) then

       call FILEIO_read( LST(:,:),                                        & ![OUT]
                         CPL_RESTART_IN_BASENAME, 'LST',     'XY', step=1 ) ![IN]
       call FILEIO_read( UST(:,:),                                        & ![OUT]
                         CPL_RESTART_IN_BASENAME, 'UST',     'XY', step=1 ) ![IN]
       call FILEIO_read( SST(:,:),                                        & ![OUT]
                         CPL_RESTART_IN_BASENAME, 'SST',     'XY', step=1 ) ![IN]
       call FILEIO_read( ALBW(:,:,I_SW),                                  & ![OUT]
                         CPL_RESTART_IN_BASENAME, 'ALBW_SW', 'XY', step=1 ) ![IN]
       call FILEIO_read( ALBW(:,:,I_LW),                                  & ![OUT]
                         CPL_RESTART_IN_BASENAME, 'ALBW_LW', 'XY', step=1 ) ![IN]
       call FILEIO_read( ALBG(:,:,I_SW),                                  & ![OUT]
                         CPL_RESTART_IN_BASENAME, 'ALBG_SW', 'XY', step=1 ) ![IN]
       call FILEIO_read( ALBG(:,:,I_LW),                                  & ![OUT]
                         CPL_RESTART_IN_BASENAME, 'ALBG_LW', 'XY', step=1 ) ![IN]
       call FILEIO_read( Z0W(:,:),                                        & ![OUT]
                         CPL_RESTART_IN_BASENAME, 'Z0W',     'XY', step=1 ) ![IN]
       call FILEIO_read( SkinT(:,:),                                      & ![OUT]
                         CPL_RESTART_IN_BASENAME, 'SkinT',   'XY', step=1 ) ![IN]
       call FILEIO_read( SkinW(:,:),                                      & ![OUT]
                         CPL_RESTART_IN_BASENAME, 'SkinW',   'XY', step=1 ) ![IN]
       call FILEIO_read( SnowQ(:,:),                                      & ![OUT]
                         CPL_RESTART_IN_BASENAME, 'SnowQ',   'XY', step=1 ) ![IN]
       call FILEIO_read( SnowT(:,:),                                      & ![OUT]
                         CPL_RESTART_IN_BASENAME, 'SnowT',   'XY', step=1 ) ![IN]

       call CPL_vars_fillhalo

       call CPL_vars_total
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** restart file for coupler is not specified.'

       LST  (:,:)   = UNDEF
       UST  (:,:)   = UNDEF
       SST  (:,:)   = UNDEF
       ALBW (:,:,:) = UNDEF
       ALBG (:,:,:) = UNDEF
       Z0W  (:,:)   = UNDEF
       SkinT(:,:)   = UNDEF
       SkinW(:,:)   = UNDEF
       SnowQ(:,:)   = UNDEF
       SnowT(:,:)   = UNDEF
    endif

    return
  end subroutine CPL_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Write coupler restart
  subroutine CPL_vars_restart_write
    use scale_const, only: &
       I_SW => CONST_I_SW, &
       I_LW => CONST_I_LW
    use scale_time, only: &
       TIME_gettimelabel
    use scale_fileio, only: &
       FILEIO_write
    implicit none

    character(len=15)     :: timelabel
    character(len=H_LONG) :: basename

    integer :: n
    !---------------------------------------------------------------------------

    if ( CPL_RESTART_OUT_BASENAME /= '' ) then

       call TIME_gettimelabel( timelabel )
       write(basename,'(A,A,A)') trim(CPL_RESTART_OUT_BASENAME), '_', trim(timelabel)

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Output restart file (coupler) ***'
       if( IO_L ) write(IO_FID_LOG,*) '*** filename: ', trim(basename)

       call CPL_vars_total

       call FILEIO_write( LST(:,:), basename,                                      CPL_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_LST), VAR_DESC(I_LST), VAR_UNIT(I_LST), 'XY', CPL_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( UST(:,:), basename,                                      CPL_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_UST), VAR_DESC(I_UST), VAR_UNIT(I_UST), 'XY', CPL_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( SST(:,:), basename,                                      CPL_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_SST), VAR_DESC(I_SST), VAR_UNIT(I_SST), 'XY', CPL_RESTART_OUT_DTYPE  ) ! [IN]

       call FILEIO_write( ALBW(:,:,I_SW),  basename, CPL_RESTART_OUT_TITLE,               & ! [IN]
                          VAR_NAME(I_ALBW_SW), VAR_DESC(I_ALBW_SW), VAR_UNIT(I_ALBW_SW),  & ! [IN]
                          'XY', CPL_RESTART_OUT_DTYPE, append=.true.                      ) ! [IN]
       call FILEIO_write( ALBW(:,:,I_LW),  basename, CPL_RESTART_OUT_TITLE,               & ! [IN]
                          VAR_NAME(I_ALBW_LW), VAR_DESC(I_ALBW_LW), VAR_UNIT(I_ALBW_LW),  & ! [IN]
                          'XY', CPL_RESTART_OUT_DTYPE, append=.true.                      ) ! [IN]
       call FILEIO_write( ALBG(:,:,I_SW),  basename, CPL_RESTART_OUT_TITLE,               & ! [IN]
                          VAR_NAME(I_ALBG_SW), VAR_DESC(I_ALBG_SW), VAR_UNIT(I_ALBG_SW),  & ! [IN]
                          'XY', CPL_RESTART_OUT_DTYPE, append=.true.                      ) ! [IN]
       call FILEIO_write( ALBG(:,:,I_LW),  basename, CPL_RESTART_OUT_TITLE,               & ! [IN]
                          VAR_NAME(I_ALBG_LW), VAR_DESC(I_ALBG_LW), VAR_UNIT(I_ALBG_LW),  & ! [IN]
                          'XY', CPL_RESTART_OUT_DTYPE, append=.true.                      ) ! [IN]
       call FILEIO_write( Z0W(:,:),        basename, CPL_RESTART_OUT_TITLE,               & ! [IN]
                          VAR_NAME(I_Z0W),     VAR_DESC(I_Z0W),     VAR_UNIT(I_Z0W),      & ! [IN]
                          'XY', CPL_RESTART_OUT_DTYPE, append=.true.                      ) ! [IN]
       call FILEIO_write( SkinT(:,:),      basename, CPL_RESTART_OUT_TITLE,               & ! [IN]
                          VAR_NAME(I_SkinT),   VAR_DESC(I_SkinT),   VAR_UNIT(I_SkinT),    & ! [IN]
                          'XY', CPL_RESTART_OUT_DTYPE, append=.true.                      ) ! [IN]
       call FILEIO_write( SkinW(:,:),      basename, CPL_RESTART_OUT_TITLE,               & ! [IN]
                          VAR_NAME(I_SkinW),   VAR_DESC(I_SkinW),   VAR_UNIT(I_SkinW),    & ! [IN]
                          'XY', CPL_RESTART_OUT_DTYPE, append=.true.                      ) ! [IN]
       call FILEIO_write( SnowQ(:,:),      basename, CPL_RESTART_OUT_TITLE,               & ! [IN]
                          VAR_NAME(I_SnowQ),   VAR_DESC(I_SnowQ),   VAR_UNIT(I_SnowQ),    & ! [IN]
                          'XY', CPL_RESTART_OUT_DTYPE, append=.true.                      ) ! [IN]
       call FILEIO_write( SnowT(:,:),      basename, CPL_RESTART_OUT_TITLE,               & ! [IN]
                          VAR_NAME(I_SnowT),   VAR_DESC(I_SnowT),   VAR_UNIT(I_SnowT),    & ! [IN]
                          'XY', CPL_RESTART_OUT_DTYPE, append=.true.                      ) ! [IN]

    endif

    return
  end subroutine CPL_vars_restart_write

  !-----------------------------------------------------------------------------
  !> History output set for coupler variables
  subroutine CPL_vars_history
    use scale_time, only: &
       TIME_DTSEC_CPL
    use scale_history, only: &
       HIST_in
    use scale_const, only: &
       I_SW => CONST_I_SW, &
       I_LW => CONST_I_LW
    implicit none
    !---------------------------------------------------------------------------

    if ( CPL_VARS_CHECKRANGE ) then
       call VALCHECK( LST  (:,:),      0.0_RP, 1000.0_RP, VAR_NAME(I_LST)    , __FILE__, __LINE__ )
       call VALCHECK( UST  (:,:),      0.0_RP, 1000.0_RP, VAR_NAME(I_UST)    , __FILE__, __LINE__ )
       call VALCHECK( SST  (:,:),      0.0_RP, 1000.0_RP, VAR_NAME(I_SST)    , __FILE__, __LINE__ )
       call VALCHECK( ALBW (:,:,I_SW), 0.0_RP,    2.0_RP, VAR_NAME(I_ALBW_SW), __FILE__, __LINE__ )
       call VALCHECK( ALBW (:,:,I_LW), 0.0_RP,    2.0_RP, VAR_NAME(I_ALBW_LW), __FILE__, __LINE__ )
       call VALCHECK( ALBG (:,:,I_SW), 0.0_RP,    2.0_RP, VAR_NAME(I_ALBG_SW), __FILE__, __LINE__ )
       call VALCHECK( ALBG (:,:,I_LW), 0.0_RP,    2.0_RP, VAR_NAME(I_ALBG_LW), __FILE__, __LINE__ )
       call VALCHECK( Z0W  (:,:),      0.0_RP, 1000.0_RP, VAR_NAME(I_Z0W)    , __FILE__, __LINE__ )
       call VALCHECK( SkinT(:,:),      0.0_RP, 1000.0_RP, VAR_NAME(I_SkinT)  , __FILE__, __LINE__ )
       call VALCHECK( SkinW(:,:),      0.0_RP, 1000.0_RP, VAR_NAME(I_SkinW)  , __FILE__, __LINE__ )
       call VALCHECK( SnowQ(:,:),      0.0_RP, 1000.0_RP, VAR_NAME(I_SnowQ)  , __FILE__, __LINE__ )
       call VALCHECK( SnowT(:,:),      0.0_RP, 1000.0_RP, VAR_NAME(I_SnowT)  , __FILE__, __LINE__ )
    endif

    call HIST_in( LST  (:,:),      'LST',      VAR_DESC(I_LST),      VAR_UNIT(I_LST),      TIME_DTSEC_CPL )
    call HIST_in( UST  (:,:),      'UST',      VAR_DESC(I_UST),      VAR_UNIT(I_UST),      TIME_DTSEC_CPL )
    call HIST_in( SST  (:,:),      'SST',      VAR_DESC(I_SST),      VAR_UNIT(I_SST),      TIME_DTSEC_CPL )
    call HIST_in( ALBW (:,:,I_SW), 'ALBW_SW',  VAR_DESC(I_ALBW_SW),  VAR_UNIT(I_ALBW_SW),  TIME_DTSEC_CPL )
    call HIST_in( ALBW (:,:,I_LW), 'ALBW_LW',  VAR_DESC(I_ALBW_LW),  VAR_UNIT(I_ALBW_LW),  TIME_DTSEC_CPL )
    call HIST_in( ALBG (:,:,I_SW), 'ALBG_SW',  VAR_DESC(I_ALBG_SW),  VAR_UNIT(I_ALBG_SW),  TIME_DTSEC_CPL )
    call HIST_in( ALBG (:,:,I_LW), 'ALBG_LW',  VAR_DESC(I_ALBG_LW),  VAR_UNIT(I_ALBG_LW),  TIME_DTSEC_CPL )
    call HIST_in( Z0W  (:,:),      'Z0W',      VAR_DESC(I_Z0W),      VAR_UNIT(I_Z0W),      TIME_DTSEC_CPL )
    call HIST_in( SkinT(:,:),      'SkinT',    VAR_DESC(I_SkinT),    VAR_UNIT(I_SkinT),    TIME_DTSEC_CPL )
    call HIST_in( SkinW(:,:),      'SkinW',    VAR_DESC(I_SkinW),    VAR_UNIT(I_SkinW),    TIME_DTSEC_CPL )
    call HIST_in( SnowQ(:,:),      'SnowQ',    VAR_DESC(I_SnowQ),    VAR_UNIT(I_SnowQ),    TIME_DTSEC_CPL )
    call HIST_in( SnowT(:,:),      'SnowT',    VAR_DESC(I_SnowT),    VAR_UNIT(I_SnowT),    TIME_DTSEC_CPL )

    return
  end subroutine CPL_vars_history

  !-----------------------------------------------------------------------------
  !> Budget monitor for coupler
  subroutine CPL_vars_total
    use scale_stats, only: &
       STAT_checktotal, &
       STAT_total
    use scale_const, only: &
       I_SW => CONST_I_SW, &
       I_LW => CONST_I_LW
    implicit none

    real(RP) :: total
    !---------------------------------------------------------------------------

    if ( STAT_checktotal ) then

!       call STAT_total( total, LST(:,:),        VAR_NAME(I_LST)     )
!       call STAT_total( total, UST(:,:),        VAR_NAME(I_UST)     )
!       call STAT_total( total, SST(:,:),        VAR_NAME(I_SST)     )
!       call STAT_total( total, ALBW(:,:,I_SW),  VAR_NAME(I_ALBW_SW) )
!       call STAT_total( total, ALBW(:,:,I_LW),  VAR_NAME(I_ALBW_LW) )
!       call STAT_total( total, ALBG(:,:,I_SW),  VAR_NAME(I_ALBG_SW) )
!       call STAT_total( total, ALBG(:,:,I_LW),  VAR_NAME(I_ALBG_LW) )
!       call STAT_total( total, Z0W(:,:),        VAR_NAME(I_Z0W)     )
!       call STAT_total( total, SkinT(:,:),      VAR_NAME(I_SkinT)   )
!       call STAT_total( total, SkinW(:,:),      VAR_NAME(I_SkinW)   )
!       call STAT_total( total, SnowQ(:,:),      VAR_NAME(I_SnowQ)   )
!       call STAT_total( total, SnowT(:,:),      VAR_NAME(I_SnowT)   )

    endif

    return
  end subroutine CPL_vars_total

  subroutine CPL_vars_merge
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    use scale_landuse, only: &
      frac_ocean => LANDUSE_frac_ocean, &
      frac_land  => LANDUSE_frac_land,  &
      frac_lake  => LANDUSE_frac_lake,  &
      frac_urban => LANDUSE_frac_urban, &
      frac_PFT   => LANDUSE_frac_PFT
    implicit none
    !---------------------------------------------------------------------------

    SkinT (:,:) = frac_ocean(:,:) * SST(:,:) &
                + frac_land (:,:) * LST(:,:) &
                + frac_urban(:,:) * UST(:,:)

    Atm_ZMFLX (:,:) = frac_ocean(:,:) * CPL_AtmOcn_ZMFLX (:,:) &
                    + frac_land (:,:) * CPL_AtmLnd_ZMFLX (:,:) &
                    + frac_urban(:,:) * CPL_AtmUrb_ZMFLX (:,:)

    Atm_XMFLX (:,:) = frac_ocean(:,:) * CPL_AtmOcn_XMFLX (:,:) &
                    + frac_land (:,:) * CPL_AtmLnd_XMFLX (:,:) &
                    + frac_urban(:,:) * CPL_AtmUrb_XMFLX (:,:)

    Atm_YMFLX (:,:) = frac_land (:,:) * CPL_AtmOcn_YMFLX (:,:) &
                    + frac_land (:,:) * CPL_AtmLnd_YMFLX (:,:) &
                    + frac_urban(:,:) * CPL_AtmUrb_YMFLX (:,:)

    Atm_SHFLX (:,:) = frac_ocean(:,:) * CPL_AtmOcn_SHFLX (:,:) &
                    + frac_land (:,:) * CPL_AtmLnd_SHFLX (:,:) &
                    + frac_urban(:,:) * CPL_AtmUrb_SHFLX (:,:)

    Atm_LHFLX (:,:) = frac_ocean(:,:) * CPL_AtmOcn_LHFLX (:,:) &
                    + frac_land (:,:) * CPL_AtmLnd_LHFLX (:,:) &
                    + frac_urban(:,:) * CPL_AtmUrb_LHFLX (:,:)

    Atm_QVFLX (:,:) = frac_ocean(:,:) * CPL_AtmOcn_QVFLX (:,:) &
                    + frac_land (:,:) * CPL_AtmLnd_QVFLX (:,:) &
                    + frac_urban(:,:) * CPL_AtmUrb_QVFLX (:,:)

    call COMM_vars8( Atm_ZMFLX(:,:), 1 )
    call COMM_vars8( Atm_XMFLX(:,:), 2 )
    call COMM_vars8( Atm_YMFLX(:,:), 3 )
    call COMM_vars8( Atm_SHFLX(:,:), 4 )
    call COMM_vars8( Atm_LHFLX(:,:), 5 )
    call COMM_vars8( Atm_QVFLX(:,:), 6 )

    call COMM_wait ( Atm_ZMFLX(:,:), 1 )
    call COMM_wait ( Atm_XMFLX(:,:), 2 )
    call COMM_wait ( Atm_YMFLX(:,:), 3 )
    call COMM_wait ( Atm_SHFLX(:,:), 4 )
    call COMM_wait ( Atm_LHFLX(:,:), 5 )
    call COMM_wait ( Atm_QVFLX(:,:), 6 )

    return
  end subroutine CPL_vars_merge

  subroutine CPL_putAtm( &
      pATM_TEMP,   &
      pATM_PRES,   &
      pATM_W,      &
      pATM_U,      &
      pATM_V,      &
      pATM_DENS,   &
      pATM_QTRC,   &
      pSFC_PRES,   &
      pSFLX_LW_dn, &
      pSFLX_SW_dn, &
      pSFLX_rain,  &
      pSFLX_snow   )
    implicit none

    real(RP), intent(in) :: pATM_TEMP  (IA,JA)
    real(RP), intent(in) :: pATM_PRES  (IA,JA)
    real(RP), intent(in) :: pATM_W     (IA,JA)
    real(RP), intent(in) :: pATM_U     (IA,JA)
    real(RP), intent(in) :: pATM_V     (IA,JA)
    real(RP), intent(in) :: pATM_DENS  (IA,JA)
    real(RP), intent(in) :: pATM_QTRC  (IA,JA,QA)
    real(RP), intent(in) :: pSFC_PRES  (IA,JA)
    real(RP), intent(in) :: pSFLX_LW_dn(IA,JA)
    real(RP), intent(in) :: pSFLX_SW_dn(IA,JA)
    real(RP), intent(in) :: pSFLX_rain (IA,JA)
    real(RP), intent(in) :: pSFLX_snow (IA,JA)
    !---------------------------------------------------------------------------

    CPL_TMPA(:,:) = pATM_TEMP  (:,:)
    CPL_PRSA(:,:) = pATM_PRES  (:,:)
    CPL_MOMZ(:,:) = pATM_W     (:,:)
    CPL_MOMX(:,:) = pATM_U     (:,:)
    CPL_MOMY(:,:) = pATM_V     (:,:)
    CPL_DENS(:,:) = pATM_DENS  (:,:)
    CPL_QV  (:,:) = pATM_QTRC  (:,:,1)
    CPL_PRSS(:,:) = pSFC_PRES  (:,:)
    CPL_LWD (:,:) = pSFLX_LW_dn(:,:)
    CPL_SWD (:,:) = pSFLX_SW_dn(:,:)
    CPL_PREC(:,:) = pSFLX_rain (:,:) &
                  + pSFLX_snow (:,:)

    return
  end subroutine CPL_putAtm

  subroutine CPL_putLnd( &
      pTG,   & ! (in)
      pQVEF, & ! (in)
      pTCS,  & ! (in)
      pDZG,  & ! (in)
      pZ0M,  & ! (in)
      pZ0H,  & ! (in)
      pZ0E   ) ! (in)
    implicit none

    real(RP), intent(in) :: pTG  (IA,JA)
    real(RP), intent(in) :: pQVEF(IA,JA)
    real(RP), intent(in) :: pTCS (IA,JA)
    real(RP), intent(in) :: pDZG (IA,JA)
    real(RP), intent(in) :: pZ0M (IA,JA)
    real(RP), intent(in) :: pZ0H (IA,JA)
    real(RP), intent(in) :: pZ0E (IA,JA)
    !---------------------------------------------------------------------------

    CPL_TG  (:,:) = pTG  (:,:)
    CPL_QVEF(:,:) = pQVEF(:,:)
    CPL_TCS (:,:) = pTCS (:,:)
    CPL_DZG (:,:) = pDZG (:,:)
    CPL_Z0M (:,:) = pZ0M (:,:)
    CPL_Z0H (:,:) = pZ0H (:,:)
    CPL_Z0E (:,:) = pZ0E (:,:)

    return
  end subroutine CPL_putLnd

  ! tentative
  subroutine CPL_putUrb
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine CPL_putUrb

  subroutine CPL_putOcn( &
      pTW ) ! (in)
    implicit none

    real(RP), intent(in) :: pTW  (IA,JA)
    !---------------------------------------------------------------------------

    CPL_TW  (:,:) = pTW  (:,:)

    return
  end subroutine CPL_putOcn

  subroutine CPL_getAtm( &
      pSFC_Z0,     &
      pSFLX_MW,    &
      pSFLX_MU,    &
      pSFLX_MV,    &
      pSFLX_SH,    &
      pSFLX_LH,    &
      pSFLX_QTRC,  &
      pUabs10,     &
      pU10,        &
      pV10,        &
      pT2,         &
      pQ2          )
    implicit none

    real(RP), intent(out) :: pSFC_Z0    (IA,JA)
    real(RP), intent(out) :: pSFLX_MW   (IA,JA)
    real(RP), intent(out) :: pSFLX_MU   (IA,JA)
    real(RP), intent(out) :: pSFLX_MV   (IA,JA)
    real(RP), intent(out) :: pSFLX_SH   (IA,JA)
    real(RP), intent(out) :: pSFLX_LH   (IA,JA)
    real(RP), intent(out) :: pSFLX_QTRC (IA,JA,QA)
    real(RP), intent(out) :: pUabs10    (IA,JA)
    real(RP), intent(out) :: pU10       (IA,JA)
    real(RP), intent(out) :: pV10       (IA,JA)
    real(RP), intent(out) :: pT2        (IA,JA)
    real(RP), intent(out) :: pQ2        (IA,JA)
    !---------------------------------------------------------------------------

    pSFC_Z0   (:,:)   = 0.0_RP ! tentative
    pSFLX_MW  (:,:)   = Atm_ZMFLX (:,:)
    pSFLX_MU  (:,:)   = Atm_XMFLX (:,:)
    pSFLX_MV  (:,:)   = Atm_YMFLX (:,:)
    pSFLX_SH  (:,:)   = Atm_SHFLX (:,:)
    pSFLX_LH  (:,:)   = Atm_LHFLX (:,:)
    pSFLX_QTRC(:,:,:) = 0.0_RP ! tentative
    pSFLX_QTRC(:,:,1) = Atm_QVFLX (:,:) ! tentative
    pUabs10   (:,:)   = 0.0_RP ! tentative
    pU10      (:,:)   = 0.0_RP ! tentative
    pV10      (:,:)   = 0.0_RP ! tentative
    pT2       (:,:)   = 0.0_RP ! tentative
    pQ2       (:,:)   = 0.0_RP ! tentative

    CNT_Atm_Lnd = 0.0_RP
    CNT_Atm_Urb = 0.0_RP
    CNT_Atm_Ocn = 0.0_RP

    return
  end subroutine CPL_getAtm

  subroutine CPL_getAtm_RD( &
      pSFC_TEMP,       &
      pSFC_albedo_land )
    implicit none

    real(RP), intent(out) :: pSFC_TEMP       (IA,JA)
    real(RP), intent(out) :: pSFC_albedo_land(IA,JA,2)
    !---------------------------------------------------------------------------

    pSFC_TEMP       (:,:)   = SkinT(:,:)
    pSFC_albedo_land(:,:,:) = ALBG (:,:,:)

    return
  end subroutine CPL_getAtm_RD

  subroutine CPL_getLnd( &
      pLnd_GHFLX,   & ! (out)
      pLnd_PRECFLX, & ! (out)
      pLnd_QVFLX    ) ! (out)
    implicit none

    real(RP), intent(out) :: pLnd_GHFLX  (IA,JA)
    real(RP), intent(out) :: pLnd_PRECFLX(IA,JA)
    real(RP), intent(out) :: pLnd_QVFLX  (IA,JA)

    pLnd_GHFLX  (:,:) = Lnd_GHFLX  (:,:)
    pLnd_PRECFLX(:,:) = Lnd_PRECFLX(:,:)
    pLnd_QVFLX  (:,:) = Lnd_QVFLX  (:,:)

    CNT_Lnd = 0.0_RP

    return
  end subroutine CPL_getLnd

  subroutine CPL_getUrb( &
      pUrb_GHFLX,   & ! (out)
      pUrb_PRECFLX, & ! (out)
      pUrb_QVFLX    ) ! (out)
    implicit none

    real(RP), intent(out) :: pUrb_GHFLX  (IA,JA)
    real(RP), intent(out) :: pUrb_PRECFLX(IA,JA)
    real(RP), intent(out) :: pUrb_QVFLX  (IA,JA)

    pUrb_GHFLX  (:,:) = Urb_GHFLX  (:,:)
    pUrb_PRECFLX(:,:) = Urb_PRECFLX(:,:)
    pUrb_QVFLX  (:,:) = Urb_QVFLX  (:,:)

    CNT_Urb = 0.0_RP

    return
  end subroutine CPL_getUrb

  subroutine CPL_getOcn( &
      pOcn_WHFLX,   & ! (out)
      pOcn_PRECFLX, & ! (out)
      pOcn_QVFLX    ) ! (out)
    implicit none

    real(RP), intent(out) :: pOcn_WHFLX  (IA,JA)
    real(RP), intent(out) :: pOcn_PRECFLX(IA,JA)
    real(RP), intent(out) :: pOcn_QVFLX  (IA,JA)

    pOcn_WHFLX  (:,:) = Ocn_WHFLX  (:,:)
    pOcn_PRECFLX(:,:) = Ocn_PRECFLX(:,:)
    pOcn_QVFLX  (:,:) = Ocn_QVFLX  (:,:)

    CNT_Ocn = 0.0_RP

    return
  end subroutine CPL_getOcn

end module mod_CPL_vars
