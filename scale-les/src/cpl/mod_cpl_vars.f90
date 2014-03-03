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
  use mod_precision
  use mod_stdio
  use mod_prof
  use mod_debug
  use mod_grid_index
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
  public :: CPL_putOcn
  public :: CPL_getCPL2Atm
  public :: CPL_getCPL2Lnd
  public :: CPL_getCPL2Ocn

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  character(len=H_SHORT), public, save :: CPL_TYPE_AtmLnd = 'OFF'   !< atmos-land coupler type
  character(len=H_SHORT), public, save :: CPL_TYPE_AtmOcn = 'OFF'   !< atmos-ocean coupler type
  logical,                public, save :: CPL_LST_UPDATE  = .false. !< is Land Surface Temperature updated?
  logical,                public, save :: CPL_SST_UPDATE  = .false. !< is Sea Surface Temperature updated?
  logical,                public, save :: CPL_sw_AtmLnd             !< do atmos-land coupler calculation?
  logical,                public, save :: CPL_sw_AtmOcn             !< do atmos-ocean coupler calculation?
  logical,                public, save :: CPL_sw_restart            !< output coupler restart?

  real(RP), public, save, allocatable :: LST  (:,:) ! Land Surface Temperature [K]
  real(RP), public, save, allocatable :: SST  (:,:) ! Sea Surface Temperature [K]
  real(RP), public, save, allocatable :: ALBW (:,:) ! Sea Surface albedo [0-1]
  real(RP), public, save, allocatable :: Z0W  (:,:) ! Sea Surface roughness length [m]
  real(RP), public, save, allocatable :: SkinT(:,:) ! Ground Skin Temperature [K]
  real(RP), public, save, allocatable :: SkinW(:,:) ! Ground Skin Water       [kg/m2]
  real(RP), public, save, allocatable :: SnowQ(:,:) ! Ground Snow amount      [kg/m2]
  real(RP), public, save, allocatable :: SnowT(:,:) ! Ground Snow Temperature [K]

  real(RP), public, save, allocatable :: XMFLX (:,:) ! x-momentum flux [kg/m2/s]
  real(RP), public, save, allocatable :: YMFLX (:,:) ! y-momentum flux [kg/m2/s]
  real(RP), public, save, allocatable :: ZMFLX (:,:) ! z-momentum flux [kg/m2/s]
  real(RP), public, save, allocatable :: SWUFLX(:,:) ! upward short-wave radiation flux (upward positive) [W/m2]
  real(RP), public, save, allocatable :: LWUFLX(:,:) ! upward long-wave radiation flux (upward positive) [W/m2]
  real(RP), public, save, allocatable :: SHFLX (:,:) ! sensible heat flux (upward positive) [W/m2]
  real(RP), public, save, allocatable :: LHFLX (:,:) ! latent heat flux (upward positive) [W/m2]
  real(RP), public, save, allocatable :: QVFLX (:,:) ! moisture flux for atmosphere [kg/m2/s]

  real(RP), public, save, allocatable :: Lnd_GHFLX  (:,:) ! ground heat flux (upward positive) [W/m2]
  real(RP), public, save, allocatable :: Lnd_PRECFLX(:,:) ! precipitation flux [kg/m2/s]
  real(RP), public, save, allocatable :: Lnd_QVFLX  (:,:) ! moisture flux for land [kg/m2/s]

  real(RP), public, save, allocatable :: Ocn_WHFLX  (:,:) ! water heat flux (upward positive) [W/m2]
  real(RP), public, save, allocatable :: Ocn_PRECFLX(:,:) ! precipitation flux [kg/m2/s]
  real(RP), public, save, allocatable :: Ocn_QVFLX  (:,:) ! moisture flux for ocean [kg/m2/s]

  ! surface fluxes from atmosphere-land coupler
  real(RP), public, save, allocatable :: AtmLnd_XMFLX      (:,:) ! x-momentum flux [kg/m2/s]
  real(RP), public, save, allocatable :: AtmLnd_YMFLX      (:,:) ! y-momentum flux [kg/m2/s]
  real(RP), public, save, allocatable :: AtmLnd_ZMFLX      (:,:) ! z-momentum flux [kg/m2/s]
  real(RP), public, save, allocatable :: AtmLnd_SWUFLX     (:,:) ! upward short-wave radiation flux (upward positive) [W/m2]
  real(RP), public, save, allocatable :: AtmLnd_LWUFLX     (:,:) ! upward long-wave radiation flux (upward positive) [W/m2]
  real(RP), public, save, allocatable :: AtmLnd_SHFLX      (:,:) ! sensible heat flux (upward positive) [W/m2]
  real(RP), public, save, allocatable :: AtmLnd_LHFLX      (:,:) ! latent heat flux (upward positive) [W/m2]
  real(RP), public, save, allocatable :: AtmLnd_QVFLX      (:,:) ! moisture flux for atmosphere [kg/m2/s]

  ! surface fluxes from atmosphere-ocean coupler
  real(RP), public, save, allocatable :: AtmOcn_XMFLX      (:,:) ! x-momentum flux [kg/m2/s]
  real(RP), public, save, allocatable :: AtmOcn_YMFLX      (:,:) ! y-momentum flux [kg/m2/s]
  real(RP), public, save, allocatable :: AtmOcn_ZMFLX      (:,:) ! z-momentum flux [kg/m2/s]
  real(RP), public, save, allocatable :: AtmOcn_SWUFLX     (:,:) ! upward short-wave radiation flux (upward positive) [W/m2]
  real(RP), public, save, allocatable :: AtmOcn_LWUFLX     (:,:) ! upward long-wave radiation flux (upward positive) [W/m2]
  real(RP), public, save, allocatable :: AtmOcn_SHFLX      (:,:) ! sensible heat flux (upward positive) [W/m2]
  real(RP), public, save, allocatable :: AtmOcn_LHFLX      (:,:) ! latent heat flux (upward positive) [W/m2]
  real(RP), public, save, allocatable :: AtmOcn_QVFLX      (:,:) ! moisture flux for atmosphere [kg/m2/s]

  ! Atmospheric values
  real(RP), public, save, allocatable :: CPL_DENS(:,:) ! air density [kg/m3]
  real(RP), public, save, allocatable :: CPL_MOMX(:,:) ! momentum x [kg/m2/s]
  real(RP), public, save, allocatable :: CPL_MOMY(:,:) ! momentum y [kg/m2/s]
  real(RP), public, save, allocatable :: CPL_MOMZ(:,:) ! momentum z [kg/m2/s]
  real(RP), public, save, allocatable :: CPL_RHOS(:,:) ! air density at the surface [kg/m3]
  real(RP), public, save, allocatable :: CPL_PRES(:,:) ! pressure at the surface [Pa]
  real(RP), public, save, allocatable :: CPL_ATMP(:,:) ! air temperature at the surface [K]
  real(RP), public, save, allocatable :: CPL_QV  (:,:) ! ratio of mass of tracer to total mass [kg/kg]
  real(RP), public, save, allocatable :: CPL_PREC(:,:) ! surface precipitation rate [kg/m2/s]
  real(RP), public, save, allocatable :: CPL_SWD (:,:) ! downward short-wave radiation flux (upward positive) [W/m2]
  real(RP), public, save, allocatable :: CPL_LWD (:,:) ! downward long-wave radiation flux (upward positive) [W/m2]

  ! Land values
  real(RP), public, save, allocatable :: CPL_TG  (:,:) ! soil temperature [K]
  real(RP), public, save, allocatable :: CPL_QVEF(:,:) ! efficiency of evaporation [no unit]
  real(RP), public, save, allocatable :: CPL_EMIT(:,:) ! emissivity in long-wave radiation [no unit]
  real(RP), public, save, allocatable :: CPL_ALBG(:,:) ! surface albedo in short-wave radiation for soil [no unit]
  real(RP), public, save, allocatable :: CPL_TCS (:,:) ! thermal conductivity for soil [W/m/K]
  real(RP), public, save, allocatable :: CPL_DZG (:,:) ! soil depth [m]
  real(RP), public, save, allocatable :: CPL_Z0M (:,:) ! roughness length for momemtum [m]
  real(RP), public, save, allocatable :: CPL_Z0H (:,:) ! roughness length for heat [m]
  real(RP), public, save, allocatable :: CPL_Z0E (:,:) ! roughness length for vapor [m]

  ! Ocean values
  real(RP), public, save, allocatable :: CPL_TW  (:,:) ! water temperature [K]

  ! counter
  real(RP), public, save :: CNT_Atm_Lnd ! counter for atmos flux by land
  real(RP), public, save :: CNT_Atm_Ocn ! counter for atmos flux by ocean
  real(RP), public, save :: CNT_Lnd     ! counter for land flux
  real(RP), public, save :: CNT_Ocn     ! counter for ocean flux

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical,                private, save :: CPL_RESTART_OUTPUT       = .false.               !< output restart file?
  character(len=H_LONG) , private, save :: CPL_RESTART_IN_BASENAME  = ''                    !< basename of the input file
  character(len=H_LONG) , private, save :: CPL_RESTART_OUT_BASENAME = 'restart_out'         !< basename of the output file
  character(len=H_MID)  , private, save :: CPL_RESTART_OUT_TITLE    = 'SCALE-LES CPL VARS.' !< title    of the output file
  character(len=H_SHORT), private, save :: CPL_RESTART_OUT_DTYPE    = 'DEFAULT'             !< REAL4 or REAL8

  logical,                private, save :: CPL_VARS_CHECKRANGE      = .false.

  integer,                private, save :: I_LST         = 1
  integer,                private, save :: I_SST         = 2
  integer,                private, save :: I_ALBW        = 3
  integer,                private, save :: I_Z0W         = 4
  integer,                private, save :: I_SkinT       = 5
  integer,                private, save :: I_SkinW       = 6
  integer,                private, save :: I_SnowQ       = 7
  integer,                private, save :: I_SnowT       = 8
  integer,                private, save :: I_XMFLX       = 9
  integer,                private, save :: I_YMFLX       = 10
  integer,                private, save :: I_ZMFLX       = 11
  integer,                private, save :: I_SWUFLX      = 12
  integer,                private, save :: I_LWUFLX      = 13
  integer,                private, save :: I_SHFLX       = 14
  integer,                private, save :: I_LHFLX       = 15
  integer,                private, save :: I_QVFLX       = 16
  integer,                private, save :: I_Lnd_GHFLX   = 17
  integer,                private, save :: I_Lnd_PRECFLX = 18
  integer,                private, save :: I_Lnd_QVFLX   = 19
  integer,                private, save :: I_Ocn_WHFLX   = 20
  integer,                private, save :: I_Ocn_PRECFLX = 21
  integer,                private, save :: I_Ocn_QVFLX   = 22

  character(len=H_SHORT), private, save :: LP_NAME(22) !< name  of the coupler variables
  character(len=H_MID),   private, save :: LP_DESC(22) !< desc. of the coupler variables
  character(len=H_SHORT), private, save :: LP_UNIT(22) !< unit  of the coupler variables

  data LP_NAME / 'LST',         &
                 'SST',         &
                 'ALBW',        &
                 'Z0W',         &
                 'SkinT',       &
                 'SkinW',       &
                 'SnowQ',       &
                 'SnowT',       &
                 'XMFLX',       &
                 'YMFLX',       &
                 'ZMFLX',       &
                 'SWUFLX',      &
                 'LWUFLX',      &
                 'SHFLX',       &
                 'LHFLX',       &
                 'QVFLX',       &
                 'Lnd_GHFLX',   &
                 'Lnd_PRECFLX', &
                 'Lnd_QVFLX',   &
                 'Ocn_WHFLX',   &
                 'Ocn_PRECFLX', &
                 'Ocn_QVFLX'    /

  data LP_DESC / 'land surface temp.',               &
                 'sea surface temp.',                &
                 'sea surface albedo',               &
                 'sea surface roughness length',     &
                 'ground skin temp.',                &
                 'ground skin water',                &
                 'ground snow amount',               &
                 'ground snow temp.',                &
                 'x-momentum flux',                  &
                 'y-momentum flux',                  &
                 'z-momentum flux',                  &
                 'upward short-wave radiation flux', &
                 'upward long-wave radiation flux',  &
                 'sensible heat flux',               &
                 'latent heat flux',                 &
                 'moisture flux for atmosphere',     &
                 'ground heat flux',                 &
                 'precipitation flux for land',      &
                 'moisture flux for land',           &
                 'water heat flux',                  &
                 'precipitation flux for ocean',     &
                 'moisture flux for ocean'           /

  data LP_UNIT / 'K',       &
                 'K',       &
                 '0-1',     &
                 'm',       &
                 'K',       &
                 'kg/m2',   &
                 'kg/m2',   &
                 'K',       &
                 'kg/m2/s', &
                 'kg/m2/s', &
                 'kg/m2/s', &
                 'W/m2',    &
                 'W/m2',    &
                 'W/m2',    &
                 'W/m2',    &
                 'kg/m2/s', &
                 'W/m2',    &
                 'kg/m2/s', &
                 'kg/m2/s', &
                 'W/m2',    &
                 'kg/m2/s', &
                 'kg/m2/s'  /

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine CPL_vars_setup
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    NAMELIST / PARAM_CPL / &
       CPL_TYPE_AtmLnd, &
       CPL_TYPE_AtmOcn, &
       CPL_LST_UPDATE,  &
       CPL_SST_UPDATE

    NAMELIST / PARAM_CPL_VARS /  &
       CPL_RESTART_IN_BASENAME,  &
       CPL_RESTART_OUTPUT,       &
       CPL_RESTART_OUT_BASENAME, &
       CPL_RESTART_OUT_TITLE,    &
       CPL_RESTART_OUT_DTYPE,    &
       CPL_VARS_CHECKRANGE

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[CPL VARS]/Categ[CPL]'


    CNT_Atm_Lnd = 0.0_RP
    CNT_Atm_Ocn = 0.0_RP
    CNT_Lnd     = 0.0_RP
    CNT_Ocn     = 0.0_RP

    allocate( LST  (IA,JA) )
    allocate( SST  (IA,JA) )
    allocate( ALBW (IA,JA) )
    allocate( Z0W  (IA,JA) )
    allocate( SkinT(IA,JA) )
    allocate( SkinW(IA,JA) )
    allocate( SnowQ(IA,JA) )
    allocate( SnowT(IA,JA) )

    allocate( XMFLX (IA,JA) )
    allocate( YMFLX (IA,JA) )
    allocate( ZMFLX (IA,JA) )
    allocate( SWUFLX(IA,JA) )
    allocate( LWUFLX(IA,JA) )
    allocate( SHFLX (IA,JA) )
    allocate( LHFLX (IA,JA) )
    allocate( QVFLX (IA,JA) )

    allocate( Lnd_GHFLX  (IA,JA) )
    allocate( Lnd_PRECFLX(IA,JA) )
    allocate( Lnd_QVFLX  (IA,JA) )

    allocate( Ocn_WHFLX  (IA,JA) )
    allocate( Ocn_PRECFLX(IA,JA) )
    allocate( Ocn_QVFLX  (IA,JA) )

    allocate( CPL_DENS(IA,JA) )
    allocate( CPL_MOMX(IA,JA) )
    allocate( CPL_MOMY(IA,JA) )
    allocate( CPL_MOMZ(IA,JA) )
    allocate( CPL_RHOS(IA,JA) )
    allocate( CPL_PRES(IA,JA) )
    allocate( CPL_ATMP(IA,JA) )
    allocate( CPL_QV  (IA,JA) )
    allocate( CPL_PREC(IA,JA) )
    allocate( CPL_SWD (IA,JA) )
    allocate( CPL_LWD (IA,JA) )

    allocate( CPL_TG  (IA,JA) )
    allocate( CPL_QVEF(IA,JA) )
    allocate( CPL_EMIT(IA,JA) )
    allocate( CPL_ALBG(IA,JA) )
    allocate( CPL_TCS (IA,JA) )
    allocate( CPL_DZG (IA,JA) )
    allocate( CPL_Z0M (IA,JA) )
    allocate( CPL_Z0H (IA,JA) )
    allocate( CPL_Z0E (IA,JA) )

    allocate( CPL_TW  (IA,JA) )

    allocate( AtmLnd_XMFLX      (IA,JA) )
    allocate( AtmLnd_YMFLX      (IA,JA) )
    allocate( AtmLnd_ZMFLX      (IA,JA) )
    allocate( AtmLnd_SWUFLX     (IA,JA) )
    allocate( AtmLnd_LWUFLX     (IA,JA) )
    allocate( AtmLnd_SHFLX      (IA,JA) )
    allocate( AtmLnd_LHFLX      (IA,JA) )
    allocate( AtmLnd_QVFLX      (IA,JA) )

    allocate( AtmOcn_XMFLX      (IA,JA) )
    allocate( AtmOcn_YMFLX      (IA,JA) )
    allocate( AtmOcn_ZMFLX      (IA,JA) )
    allocate( AtmOcn_SWUFLX     (IA,JA) )
    allocate( AtmOcn_LWUFLX     (IA,JA) )
    allocate( AtmOcn_SHFLX      (IA,JA) )
    allocate( AtmOcn_LHFLX      (IA,JA) )
    allocate( AtmOcn_QVFLX      (IA,JA) )


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

    if( IO_L ) write(IO_FID_LOG,*) '*** [CPL] selected components'

    ! Atoms-Land Switch
    if ( CPL_TYPE_AtmLnd /= 'OFF' .AND. CPL_TYPE_AtmLnd /= 'NONE' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Atmos-Land Coupler : ON'
       CPL_sw_AtmLnd = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Atmos-Land Coupler : OFF'
       CPL_sw_AtmLnd = .false.
    endif

    ! Atoms-Ocean Switch
    if ( CPL_TYPE_AtmOcn /= 'OFF' .AND. CPL_TYPE_AtmOcn /= 'NONE' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Atmos-Ocean Coupler : ON'
       CPL_sw_AtmOcn = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Atmos-Ocean Coupler : OFF'
       CPL_sw_AtmOcn = .false.
    endif

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[CPL VARS]/Categ[CPL]'

    return
  end subroutine CPL_vars_setup

  !-----------------------------------------------------------------------------
  !> fill HALO region of coupler variables
  subroutine CPL_vars_fillhalo
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none
    !---------------------------------------------------------------------------

    ! fill IHALO & JHALO
    call COMM_vars8( LST  (:,:), 1 )
    call COMM_vars8( SST  (:,:), 2 )
    call COMM_vars8( ALBW (:,:), 3 )
    call COMM_vars8( Z0W  (:,:), 4 )
    call COMM_vars8( SkinT(:,:), 5 )
    call COMM_vars8( SkinW(:,:), 6 )
    call COMM_vars8( SnowQ(:,:), 7 )
    call COMM_vars8( SnowT(:,:), 8 )

    call COMM_wait ( LST  (:,:), 1 )
    call COMM_wait ( SST  (:,:), 2 )
    call COMM_wait ( ALBW (:,:), 3 )
    call COMM_wait ( Z0W  (:,:), 4 )
    call COMM_wait ( SkinT(:,:), 5 )
    call COMM_wait ( SkinW(:,:), 6 )
    call COMM_wait ( SnowQ(:,:), 7 )
    call COMM_wait ( SnowT(:,:), 8 )

    return
  end subroutine CPL_vars_fillhalo

  !-----------------------------------------------------------------------------
  !> Read coupler restart
  subroutine CPL_vars_restart_read
    use mod_fileio, only: &
       FILEIO_read
    use mod_const, only: &
       CONST_UNDEF
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input restart file (coupler) ***'

    call PROF_rapstart('FILE I NetCDF')

    if ( CPL_RESTART_IN_BASENAME /= '' ) then

       call FILEIO_read( LST(:,:),                                      & ![OUT]
                         CPL_RESTART_IN_BASENAME, 'LST',   'XY', step=1 ) ![IN]
       call FILEIO_read( SST(:,:),                                      & ![OUT]
                         CPL_RESTART_IN_BASENAME, 'SST',   'XY', step=1 ) ![IN]
       call FILEIO_read( ALBW(:,:),                                     & ![OUT]
                         CPL_RESTART_IN_BASENAME, 'ALBW',  'XY', step=1 ) ![IN]
       call FILEIO_read( Z0W(:,:),                                      & ![OUT]
                         CPL_RESTART_IN_BASENAME, 'Z0W',   'XY', step=1 ) ![IN]
       call FILEIO_read( SkinT(:,:),                                    & ![OUT]
                         CPL_RESTART_IN_BASENAME, 'SkinT', 'XY', step=1 ) ![IN]
       call FILEIO_read( SkinW(:,:),                                    & ![OUT]
                         CPL_RESTART_IN_BASENAME, 'SkinW', 'XY', step=1 ) ![IN]
       call FILEIO_read( SnowQ(:,:),                                    & ![OUT]
                         CPL_RESTART_IN_BASENAME, 'SnowQ', 'XY', step=1 ) ![IN]
       call FILEIO_read( SnowT(:,:),                                    & ![OUT]
                         CPL_RESTART_IN_BASENAME, 'SnowT', 'XY', step=1 ) ![IN]

       call CPL_vars_fillhalo

       call CPL_vars_total

    else

       if( IO_L ) write(IO_FID_LOG,*) '*** restart file for coupler is not specified.'

       LST  (:,:) = 300.0_RP
       SST  (:,:) = 300.0_RP
       ALBW (:,:) = 0.05_RP
       Z0W  (:,:) = 1.0E-4_RP
       SkinT(:,:) = 300.0_RP
!       LST  (:,:) = CONST_UNDEF
!       SST  (:,:) = CONST_UNDEF
!       ALBW (:,:) = CONST_UNDEF
!       Z0W  (:,:) = CONST_UNDEF
!       SkinT(:,:) = CONST_UNDEF
       SkinW(:,:) = CONST_UNDEF
       SnowQ(:,:) = CONST_UNDEF
       SnowT(:,:) = CONST_UNDEF

    endif

    call PROF_rapend  ('FILE I NetCDF')

    return
  end subroutine CPL_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Write coupler restart
  subroutine CPL_vars_restart_write
    use mod_time, only: &
       NOWSEC => TIME_NOWDAYSEC
    use mod_fileio, only: &
       FILEIO_write
    implicit none

    character(len=H_LONG) :: bname

    integer :: n
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE O NetCDF')

    if ( CPL_RESTART_OUT_BASENAME /= '' ) then

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Output restart file (coupler) ***'

       bname = ''
       write(bname(1:15), '(F15.3)') NOWSEC
       do n = 1, 15
          if ( bname(n:n) == ' ' ) bname(n:n) = '0'
       enddo
       write(bname,'(A,A,A)') trim(CPL_RESTART_OUT_BASENAME), '_', trim(bname)

       call FILEIO_write( LST(:,:),   bname, CPL_RESTART_OUT_TITLE,                                          & ! [IN]
                          LP_NAME(I_LST),   LP_DESC(I_LST),   LP_UNIT(I_LST),   'XY', CPL_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( SST(:,:),   bname, CPL_RESTART_OUT_TITLE,                                          & ! [IN]
                          LP_NAME(I_SST),   LP_DESC(I_SST),   LP_UNIT(I_SST),   'XY', CPL_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( ALBW(:,:),  bname, CPL_RESTART_OUT_TITLE,                                          & ! [IN]
                          LP_NAME(I_ALBW),  LP_DESC(I_ALBW),  LP_UNIT(I_ALBW),  'XY', CPL_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( Z0W(:,:),   bname, CPL_RESTART_OUT_TITLE,                                          & ! [IN]
                          LP_NAME(I_Z0W),   LP_DESC(I_Z0W),   LP_UNIT(I_Z0W),   'XY', CPL_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( SkinT(:,:), bname, CPL_RESTART_OUT_TITLE,                                          & ! [IN]
                          LP_NAME(I_SkinT), LP_DESC(I_SkinT), LP_UNIT(I_SkinT), 'XY', CPL_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( SkinW(:,:), bname, CPL_RESTART_OUT_TITLE,                                          & ! [IN]
                          LP_NAME(I_SkinW), LP_DESC(I_SkinW), LP_UNIT(I_SkinW), 'XY', CPL_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( SnowQ(:,:), bname, CPL_RESTART_OUT_TITLE,                                          & ! [IN]
                          LP_NAME(I_SnowQ), LP_DESC(I_SnowQ), LP_UNIT(I_SnowQ), 'XY', CPL_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( SnowT(:,:), bname, CPL_RESTART_OUT_TITLE,                                          & ! [IN]
                          LP_NAME(I_SnowT), LP_DESC(I_SnowT), LP_UNIT(I_SnowT), 'XY', CPL_RESTART_OUT_DTYPE  ) ! [IN]

    endif

    call PROF_rapend  ('FILE O NetCDF')

    call CPL_vars_total

    return
  end subroutine CPL_vars_restart_write

  !-----------------------------------------------------------------------------
  !> History output set for coupler variables
  subroutine CPL_vars_history
    use mod_time, only: &
       TIME_DTSEC_CPL
    use mod_history, only: &
       HIST_in
    implicit none
    !---------------------------------------------------------------------------

    if ( CPL_VARS_CHECKRANGE ) then
       call VALCHECK( LST  (:,:), 0.0_RP, 1000.0_RP, LP_NAME(I_LST)  , __FILE__, __LINE__ )
       call VALCHECK( SST  (:,:), 0.0_RP, 1000.0_RP, LP_NAME(I_SST)  , __FILE__, __LINE__ )
       call VALCHECK( ALBW (:,:), 0.0_RP,    2.0_RP, LP_NAME(I_ALBW) , __FILE__, __LINE__ )
       call VALCHECK( Z0W  (:,:), 0.0_RP, 1000.0_RP, LP_NAME(I_Z0W)  , __FILE__, __LINE__ )
       call VALCHECK( SkinT(:,:), 0.0_RP, 1000.0_RP, LP_NAME(I_SkinT), __FILE__, __LINE__ )
       call VALCHECK( SkinW(:,:), 0.0_RP, 1000.0_RP, LP_NAME(I_SkinW), __FILE__, __LINE__ )
       call VALCHECK( SnowQ(:,:), 0.0_RP, 1000.0_RP, LP_NAME(I_SnowQ), __FILE__, __LINE__ )
       call VALCHECK( SnowT(:,:), 0.0_RP, 1000.0_RP, LP_NAME(I_SnowT), __FILE__, __LINE__ )

       call VALCHECK( XMFLX(:,:),  -1.0E4_RP, 1.0E4_RP, LP_NAME(I_XMFLX) , __FILE__, __LINE__ )
       call VALCHECK( YMFLX(:,:),  -1.0E4_RP, 1.0E4_RP, LP_NAME(I_YMFLX) , __FILE__, __LINE__ )
       call VALCHECK( ZMFLX(:,:),  -1.0E4_RP, 1.0E4_RP, LP_NAME(I_ZMFLX) , __FILE__, __LINE__ )
       call VALCHECK( SWUFLX(:,:), -1.0E4_RP, 1.0E4_RP, LP_NAME(I_SWUFLX), __FILE__, __LINE__ )
       call VALCHECK( LWUFLX(:,:), -1.0E4_RP, 1.0E4_RP, LP_NAME(I_LWUFLX), __FILE__, __LINE__ )
       call VALCHECK( SHFLX(:,:),  -1.0E4_RP, 1.0E4_RP, LP_NAME(I_SHFLX) , __FILE__, __LINE__ )
       call VALCHECK( LHFLX(:,:),  -1.0E4_RP, 1.0E4_RP, LP_NAME(I_LHFLX) , __FILE__, __LINE__ )
       call VALCHECK( QVFLX(:,:),  -1.0E4_RP, 1.0E4_RP, LP_NAME(I_QVFLX) , __FILE__, __LINE__ )

       call VALCHECK( Lnd_GHFLX  (:,:), -1.0E4_RP, 1.0E4_RP, LP_NAME(I_Lnd_GHFLX)  , __FILE__, __LINE__ )
       call VALCHECK( Lnd_PRECFLX(:,:), -1.0E4_RP, 1.0E4_RP, LP_NAME(I_Lnd_PRECFLX), __FILE__, __LINE__ )
       call VALCHECK( Lnd_QVFLX  (:,:), -1.0E4_RP, 1.0E4_RP, LP_NAME(I_Lnd_QVFLX)  , __FILE__, __LINE__ )

       call VALCHECK( Ocn_WHFLX  (:,:), -1.0E4_RP, 1.0E4_RP, LP_NAME(I_Ocn_WHFLX)  , __FILE__, __LINE__ )
       call VALCHECK( Ocn_PRECFLX(:,:), -1.0E4_RP, 1.0E4_RP, LP_NAME(I_Ocn_PRECFLX), __FILE__, __LINE__ )
       call VALCHECK( Ocn_QVFLX  (:,:), -1.0E4_RP, 1.0E4_RP, LP_NAME(I_Ocn_QVFLX)  , __FILE__, __LINE__ )
    endif

    call HIST_in( LST  (:,:), 'LST',   LP_DESC(I_LST),   LP_UNIT(I_LST),   TIME_DTSEC_CPL )
    call HIST_in( SST  (:,:), 'SST',   LP_DESC(I_SST),   LP_UNIT(I_SST),   TIME_DTSEC_CPL )
    call HIST_in( ALBW (:,:), 'ALBW',  LP_DESC(I_ALBW),  LP_UNIT(I_ALBW),  TIME_DTSEC_CPL )
    call HIST_in( Z0W  (:,:), 'Z0W',   LP_DESC(I_Z0W),   LP_UNIT(I_Z0W),   TIME_DTSEC_CPL )
    call HIST_in( SkinT(:,:), 'SkinT', LP_DESC(I_SkinT), LP_UNIT(I_SkinT), TIME_DTSEC_CPL )
    call HIST_in( SkinW(:,:), 'SkinW', LP_DESC(I_SkinW), LP_UNIT(I_SkinW), TIME_DTSEC_CPL )
    call HIST_in( SnowQ(:,:), 'SnowQ', LP_DESC(I_SnowQ), LP_UNIT(I_SnowQ), TIME_DTSEC_CPL )
    call HIST_in( SnowT(:,:), 'SnowT', LP_DESC(I_SnowT), LP_UNIT(I_SnowT), TIME_DTSEC_CPL )

    call HIST_in( XMFLX (:,:), 'XMFLX',  LP_DESC(I_XMFLX),  LP_UNIT(I_XMFLX),  TIME_DTSEC_CPL )
    call HIST_in( YMFLX (:,:), 'YMFLX',  LP_DESC(I_YMFLX),  LP_UNIT(I_YMFLX),  TIME_DTSEC_CPL )
    call HIST_in( ZMFLX (:,:), 'ZMFLX',  LP_DESC(I_ZMFLX),  LP_UNIT(I_ZMFLX),  TIME_DTSEC_CPL )
    call HIST_in( SWUFLX(:,:), 'SWUFLX', LP_DESC(I_SWUFLX), LP_UNIT(I_SWUFLX), TIME_DTSEC_CPL )
    call HIST_in( LWUFLX(:,:), 'LWUFLX', LP_DESC(I_LWUFLX), LP_UNIT(I_LWUFLX), TIME_DTSEC_CPL )
    call HIST_in( SHFLX (:,:), 'SHFLX',  LP_DESC(I_SHFLX),  LP_UNIT(I_SHFLX),  TIME_DTSEC_CPL )
    call HIST_in( LHFLX (:,:), 'LHFLX',  LP_DESC(I_LHFLX),  LP_UNIT(I_LHFLX),  TIME_DTSEC_CPL )
    call HIST_in( QVFLX (:,:), 'QVFLX',  LP_DESC(I_QVFLX),  LP_UNIT(I_QVFLX),  TIME_DTSEC_CPL )

    call HIST_in( Lnd_GHFLX  (:,:), 'Lnd_GHFLX',   LP_DESC(I_Lnd_GHFLX),   LP_UNIT(I_Lnd_GHFLX),   TIME_DTSEC_CPL )
    call HIST_in( Lnd_PRECFLX(:,:), 'Lnd_PRECFLX', LP_DESC(I_Lnd_PRECFLX), LP_UNIT(I_Lnd_PRECFLX), TIME_DTSEC_CPL )
    call HIST_in( Lnd_QVFLX  (:,:), 'Lnd_QVFLX',   LP_DESC(I_Lnd_QVFLX),   LP_UNIT(I_Lnd_QVFLX),   TIME_DTSEC_CPL )

    call HIST_in( Ocn_WHFLX  (:,:), 'Ocn_WHFLX',   LP_DESC(I_Ocn_WHFLX),   LP_UNIT(I_Ocn_WHFLX),   TIME_DTSEC_CPL )
    call HIST_in( Ocn_PRECFLX(:,:), 'Ocn_PRECFLX', LP_DESC(I_Ocn_PRECFLX), LP_UNIT(I_Ocn_PRECFLX), TIME_DTSEC_CPL )
    call HIST_in( Ocn_QVFLX  (:,:), 'Ocn_QVFLX',   LP_DESC(I_Ocn_QVFLX),   LP_UNIT(I_Ocn_QVFLX),   TIME_DTSEC_CPL )

    return
  end subroutine CPL_vars_history

  !-----------------------------------------------------------------------------
  !> Budget monitor for coupler
  subroutine CPL_vars_total
    use mod_stats, only: &
       STAT_checktotal, &
       STAT_total
    implicit none

    real(RP) :: total
    !---------------------------------------------------------------------------

    if ( STAT_checktotal ) then

!       call STAT_total( total, LST(:,:),   LP_NAME(I_LST)   )
!       call STAT_total( total, SST(:,:),   LP_NAME(I_SST)   )
!       call STAT_total( total, ALBW(:,:),  LP_NAME(I_ALBW)  )
!       call STAT_total( total, Z0W(:,:),   LP_NAME(I_Z0W)   )
!       call STAT_total( total, SkinT(:,:), LP_NAME(I_SkinT) )
!       call STAT_total( total, SkinW(:,:), LP_NAME(I_SkinW) )
!       call STAT_total( total, SnowQ(:,:), LP_NAME(I_SnowQ) )
!       call STAT_total( total, SnowT(:,:), LP_NAME(I_SnowT) )

    endif

    return
  end subroutine CPL_vars_total

  subroutine CPL_vars_merge
    use mod_landuse, only: &
      frac_ocean => LANDUSE_frac_ocean
    implicit none
    !---------------------------------------------------------------------------

    ! tentative
    if( CPL_sw_AtmLnd ) then
      frac_ocean(:,:) = 0.0_RP
    end if

    ! input land
    SkinT (:,:) = LST          (:,:) * ( 1.0_RP - frac_ocean(:,:) )
    XMFLX (:,:) = AtmLnd_XMFLX (:,:) * ( 1.0_RP - frac_ocean(:,:) )
    YMFLX (:,:) = AtmLnd_YMFLX (:,:) * ( 1.0_RP - frac_ocean(:,:) )
    ZMFLX (:,:) = AtmLnd_ZMFLX (:,:) * ( 1.0_RP - frac_ocean(:,:) )
    SWUFLX(:,:) = AtmLnd_SWUFLX(:,:) * ( 1.0_RP - frac_ocean(:,:) )
    LWUFLX(:,:) = AtmLnd_LWUFLX(:,:) * ( 1.0_RP - frac_ocean(:,:) )
    SHFLX (:,:) = AtmLnd_SHFLX (:,:) * ( 1.0_RP - frac_ocean(:,:) )
    LHFLX (:,:) = AtmLnd_LHFLX (:,:) * ( 1.0_RP - frac_ocean(:,:) )
    QVFLX (:,:) = AtmLnd_QVFLX (:,:) * ( 1.0_RP - frac_ocean(:,:) )

    ! input ocean
    SkinT (:,:) = SkinT (:,:) + SST          (:,:) * frac_ocean(:,:)
    XMFLX (:,:) = XMFLX (:,:) + AtmOcn_XMFLX (:,:) * frac_ocean(:,:)
    YMFLX (:,:) = YMFLX (:,:) + AtmOcn_YMFLX (:,:) * frac_ocean(:,:)
    ZMFLX (:,:) = ZMFLX (:,:) + AtmOcn_ZMFLX (:,:) * frac_ocean(:,:)
    SWUFLX(:,:) = SWUFLX(:,:) + AtmOcn_SWUFLX(:,:) * frac_ocean(:,:)
    LWUFLX(:,:) = LWUFLX(:,:) + AtmOcn_LWUFLX(:,:) * frac_ocean(:,:)
    SHFLX (:,:) = SHFLX (:,:) + AtmOcn_SHFLX (:,:) * frac_ocean(:,:)
    LHFLX (:,:) = LHFLX (:,:) + AtmOcn_LHFLX (:,:) * frac_ocean(:,:)
    QVFLX (:,:) = QVFLX (:,:) + AtmOcn_QVFLX (:,:) * frac_ocean(:,:)

    return
  end subroutine CPL_vars_merge

  subroutine CPL_putAtm( &
      pDENS, pMOMX, pMOMY, pMOMZ, & ! (in)
      pRHOS, pPRES, pATMP, pQV,   & ! (in)
      pPREC, pSWD, pLWD           ) ! (in)
    implicit none

    real(RP), intent(in) :: pDENS(IA,JA)
    real(RP), intent(in) :: pMOMX(IA,JA)
    real(RP), intent(in) :: pMOMY(IA,JA)
    real(RP), intent(in) :: pMOMZ(IA,JA)
    real(RP), intent(in) :: pRHOS(IA,JA)
    real(RP), intent(in) :: pPRES(IA,JA)
    real(RP), intent(in) :: pATMP(IA,JA)
    real(RP), intent(in) :: pQV  (IA,JA)
    real(RP), intent(in) :: pPREC(IA,JA)
    real(RP), intent(in) :: pSWD (IA,JA)
    real(RP), intent(in) :: pLWD (IA,JA)
    !---------------------------------------------------------------------------

    CPL_DENS(:,:) = pDENS(:,:)
    CPL_MOMX(:,:) = pMOMX(:,:)
    CPL_MOMY(:,:) = pMOMY(:,:)
    CPL_MOMZ(:,:) = pMOMZ(:,:)
    CPL_RHOS(:,:) = pRHOS(:,:)
    CPL_PRES(:,:) = pPRES(:,:)
    CPL_ATMP(:,:) = pATMP(:,:)
    CPL_QV  (:,:) = pQV  (:,:)
    CPL_PREC(:,:) = pPREC(:,:)
    CPL_SWD (:,:) = pSWD (:,:)
    CPL_LWD (:,:) = pLWD (:,:)

    return
  end subroutine CPL_putAtm

  subroutine CPL_putLnd( &
      pTG, pQVEF, pEMIT, & ! (in)
      pALBG, pTCS, pDZG, & ! (in)
      pZ0M, pZ0H, pZ0E   ) ! (in)
    implicit none

    real(RP), intent(in) :: pTG  (IA,JA)
    real(RP), intent(in) :: pQVEF(IA,JA)
    real(RP), intent(in) :: pEMIT(IA,JA)
    real(RP), intent(in) :: pALBG(IA,JA)
    real(RP), intent(in) :: pTCS (IA,JA)
    real(RP), intent(in) :: pDZG (IA,JA)
    real(RP), intent(in) :: pZ0M (IA,JA)
    real(RP), intent(in) :: pZ0H (IA,JA)
    real(RP), intent(in) :: pZ0E (IA,JA)
    !---------------------------------------------------------------------------

    CPL_TG  (:,:) = pTG  (:,:)
    CPL_QVEF(:,:) = pQVEF(:,:)
    CPL_EMIT(:,:) = pEMIT(:,:)
    CPL_ALBG(:,:) = pALBG(:,:)
    CPL_TCS (:,:) = pTCS (:,:)
    CPL_DZG (:,:) = pDZG (:,:)
    CPL_Z0M (:,:) = pZ0M (:,:)
    CPL_Z0H (:,:) = pZ0H (:,:)
    CPL_Z0E (:,:) = pZ0E (:,:)

    return
  end subroutine CPL_putLnd

  subroutine CPL_putOcn( &
      pTW ) ! (in)
    implicit none

    real(RP), intent(in) :: pTW  (IA,JA)
    !---------------------------------------------------------------------------

    CPL_TW  (:,:) = pTW  (:,:)

    return
  end subroutine CPL_putOcn

  subroutine CPL_getCPL2Atm( &
      pXMFLX, pYMFLX, pZMFLX, & ! (out)
      pSWUFLX, pLWUFLX,       & ! (out)
      pSHFLX, pLHFLX, pQVFLX  ) ! (out)
    implicit none

    real(RP), intent(out) :: pXMFLX (IA,JA)
    real(RP), intent(out) :: pYMFLX (IA,JA)
    real(RP), intent(out) :: pZMFLX (IA,JA)
    real(RP), intent(out) :: pSWUFLX(IA,JA)
    real(RP), intent(out) :: pLWUFLX(IA,JA)
    real(RP), intent(out) :: pSHFLX (IA,JA)
    real(RP), intent(out) :: pLHFLX (IA,JA)
    real(RP), intent(out) :: pQVFLX (IA,JA)
    !---------------------------------------------------------------------------

    pXMFLX (:,:) = XMFLX (:,:)
    pYMFLX (:,:) = YMFLX (:,:)
    pZMFLX (:,:) = ZMFLX (:,:)
    pSWUFLX(:,:) = SWUFLX(:,:)
    pLWUFLX(:,:) = LWUFLX(:,:)
    pSHFLX (:,:) = SHFLX (:,:)
    pLHFLX (:,:) = LHFLX (:,:)
    pQVFLX (:,:) = QVFLX (:,:)

    CNT_Atm_Lnd = 0.0_RP
    CNT_Atm_Ocn = 0.0_RP

    return
  end subroutine CPL_getCPL2Atm

  subroutine CPL_getCPL2Lnd( &
      pLnd_GHFLX, pLnd_PRECFLX, pLnd_QVFLX ) ! (out)
    implicit none

    real(RP), intent(out) :: pLnd_GHFLX  (IA,JA)
    real(RP), intent(out) :: pLnd_PRECFLX(IA,JA)
    real(RP), intent(out) :: pLnd_QVFLX  (IA,JA)

    pLnd_GHFLX  (:,:) = Lnd_GHFLX  (:,:)
    pLnd_PRECFLX(:,:) = Lnd_PRECFLX(:,:)
    pLnd_QVFLX  (:,:) = Lnd_QVFLX  (:,:)

    CNT_Lnd = 0.0_RP

    return
  end subroutine CPL_getCPL2Lnd

  subroutine CPL_getCPL2Ocn( &
      pOcn_WHFLX, pOcn_PRECFLX, pOcn_QVFLX ) ! (out)
    implicit none

    real(RP), intent(out) :: pOcn_WHFLX  (IA,JA)
    real(RP), intent(out) :: pOcn_PRECFLX(IA,JA)
    real(RP), intent(out) :: pOcn_QVFLX  (IA,JA)

    pOcn_WHFLX  (:,:) = Ocn_WHFLX  (:,:)
    pOcn_PRECFLX(:,:) = Ocn_PRECFLX(:,:)
    pOcn_QVFLX  (:,:) = Ocn_QVFLX  (:,:)

    CNT_Ocn = 0.0_RP

    return
  end subroutine CPL_getCPL2Ocn

end module mod_CPL_vars
