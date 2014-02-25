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
  use mod_tracer
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
  public :: CPL_AtmLnd_putCPL
  public :: CPL_getCPL2Atm
  public :: CPL_getCPL2Lnd
  public :: CPL_AtmLnd_getAtm2CPL
  public :: CPL_AtmLnd_getLnd2CPL
  public :: CPL_flushAtm
  public :: CPL_flushLnd
  public :: CPL_AtmLnd_flushCPL

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, allocatable :: LST  (:,:) ! Land Surface Temperature [K]
  real(RP), public, allocatable :: SST  (:,:) ! Sea Surface Temperature [K]
  real(RP), public, allocatable :: SkinT(:,:) ! Ground Skin Temperature [K]
  real(RP), public, allocatable :: SkinW(:,:) ! Ground Skin Water       [kg/m2]
  real(RP), public, allocatable :: SnowQ(:,:) ! Ground Snow amount      [kg/m2]
  real(RP), public, allocatable :: SnowT(:,:) ! Ground Snow Temperature [K]

  character(len=H_SHORT), public, save :: CPL_TYPE_AtmLnd = 'OFF'   !< atmos-land coupler type
  character(len=H_SHORT), public, save :: CPL_TYPE_AtmOcn = 'OFF'   !< atmos-ocean coupler type
  logical,                public, save :: CPL_sw_AtmLnd             !< do atmos-land coupler calculation?
  logical,                public, save :: CPL_sw_AtmOcn             !< do atmos-ocean coupler calculation?
  logical,                public, save :: CPL_sw_restart            !< output coupler restart?
  logical,                public, save :: CPL_LST_UPDATE  = .false. !< is Land Surface Temperature updated?

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

  ! surface fluxes for atmosphere
  real(RP), private, allocatable :: XMFLX    (:,:) ! x-momentum flux [kg/m2/s]
  real(RP), private, allocatable :: YMFLX    (:,:) ! y-momentum flux [kg/m2/s]
  real(RP), private, allocatable :: ZMFLX    (:,:) ! z-momentum flux [kg/m2/s]
  real(RP), private, allocatable :: SWUFLX   (:,:) ! upward short-wave radiation flux (upward positive) [W/m2]
  real(RP), private, allocatable :: LWUFLX   (:,:) ! upward long-wave radiation flux (upward positive) [W/m2]
  real(RP), private, allocatable :: SHFLX    (:,:) ! sensible heat flux (upward positive) [W/m2]
  real(RP), private, allocatable :: LHFLX    (:,:) ! latent heat flux (upward positive) [W/m2]
  real(RP), private, allocatable :: QVFLX_Atm(:,:) ! moisture flux for atmosphere [kg/m2/s]

  ! surface fluxes for land
  real(RP), private, allocatable :: GHFLX    (:,:) ! ground heat flux (upward positive) [W/m2]
  real(RP), private, allocatable :: PRECFLX  (:,:) ! precipitation flux [kg/m2/s]
  real(RP), private, allocatable :: QVFLX_Lnd(:,:) ! moisture flux for land [kg/m2/s]

  ! Atmospheric values
  real(RP), private, allocatable :: DENS(:,:) ! air density [kg/m3]
  real(RP), private, allocatable :: MOMX(:,:) ! momentum x [kg/m2/s]
  real(RP), private, allocatable :: MOMY(:,:) ! momentum y [kg/m2/s]
  real(RP), private, allocatable :: MOMZ(:,:) ! momentum z [kg/m2/s]
  real(RP), private, allocatable :: RHOS(:,:) ! air density at the surface [kg/m3]
  real(RP), private, allocatable :: PRES(:,:) ! pressure at the surface [Pa]
  real(RP), private, allocatable :: ATMP(:,:) ! air temperature at the surface [K]
  real(RP), private, allocatable :: QV  (:,:) ! ratio of mass of tracer to total mass [kg/kg]
  real(RP), private, allocatable :: PREC(:,:) ! surface precipitation rate [kg/m2/s]
  real(RP), private, allocatable :: SWD (:,:) ! downward short-wave radiation flux (upward positive) [W/m2]
  real(RP), private, allocatable :: LWD (:,:) ! downward long-wave radiation flux (upward positive) [W/m2]

  ! Land values
  real(RP), private, allocatable :: TG  (:,:) ! soil temperature [K]
  real(RP), private, allocatable :: QVEF(:,:) ! efficiency of evaporation [no unit]
  real(RP), private, allocatable :: EMIT(:,:) ! emissivity in long-wave radiation [no unit]
  real(RP), private, allocatable :: ALB (:,:) ! surface albedo in short-wave radiation [no unit]
  real(RP), private, allocatable :: TCS (:,:) ! thermal conductivity for soil [W/m/K]
  real(RP), private, allocatable :: DZG (:,:) ! soil depth [m]
  real(RP), private, allocatable :: Z0M (:,:) ! roughness length for momemtum [m]
  real(RP), private, allocatable :: Z0H (:,:) ! roughness length for heat [m]
  real(RP), private, allocatable :: Z0E (:,:) ! roughness length for vapor [m]

  ! AtmLnd surface fluxes for atmosphere
  real(RP), private, allocatable :: Lnd_XMFLX    (:,:) ! x-momentum flux [kg/m2/s]
  real(RP), private, allocatable :: Lnd_YMFLX    (:,:) ! y-momentum flux [kg/m2/s]
  real(RP), private, allocatable :: Lnd_ZMFLX    (:,:) ! z-momentum flux [kg/m2/s]
  real(RP), private, allocatable :: Lnd_SWUFLX   (:,:) ! upward short-wave radiation flux (upward positive) [W/m2]
  real(RP), private, allocatable :: Lnd_LWUFLX   (:,:) ! upward long-wave radiation flux (upward positive) [W/m2]
  real(RP), private, allocatable :: Lnd_SHFLX    (:,:) ! sensible heat flux (upward positive) [W/m2]
  real(RP), private, allocatable :: Lnd_LHFLX    (:,:) ! latent heat flux (upward positive) [W/m2]
  real(RP), private, allocatable :: Lnd_QVFLX_Atm(:,:) ! moisture flux for atmosphere [kg/m2/s]

  ! AtmLnd surface fluxes for land
  real(RP), private, allocatable :: Lnd_GHFLX    (:,:) ! ground heat flux (upward positive) [W/m2]
  real(RP), private, allocatable :: Lnd_PRECFLX  (:,:) ! precipitation flux [kg/m2/s]
  real(RP), private, allocatable :: Lnd_QVFLX_Lnd(:,:) ! moisture flux for land [kg/m2/s]

  ! counter
  real(RP), private, save :: CNT_putAtm     ! counter for putAtm
  real(RP), private, save :: CNT_putLnd     ! counter for putLnd
  real(RP), private, save :: CNT_getCPL2Atm ! counter for getDat2Atm
  real(RP), private, save :: CNT_getCPL2Lnd ! counter for getDat2Lnd

  integer,                    private, save :: I_LST       = 1
  integer,                    private, save :: I_SST       = 2
  integer,                    private, save :: I_SkinT     = 3
  integer,                    private, save :: I_SkinW     = 4
  integer,                    private, save :: I_SnowQ     = 5
  integer,                    private, save :: I_SnowT     = 6
  integer,                    private, save :: I_XMFLX     = 7
  integer,                    private, save :: I_YMFLX     = 8
  integer,                    private, save :: I_ZMFLX     = 9
  integer,                    private, save :: I_SWUFLX    = 10
  integer,                    private, save :: I_LWUFLX    = 11
  integer,                    private, save :: I_SHFLX     = 12
  integer,                    private, save :: I_LHFLX     = 13
  integer,                    private, save :: I_QVFLX_Atm = 14
  integer,                    private, save :: I_GHFLX     = 15
  integer,                    private, save :: I_PRECFLX   = 16
  integer,                    private, save :: I_QVFLX_Lnd = 17

  character(len=H_SHORT), private, save :: LP_NAME(17) !< name  of the coupler variables
  character(len=H_MID),   private, save :: LP_DESC(17) !< desc. of the coupler variables
  character(len=H_SHORT), private, save :: LP_UNIT(17) !< unit  of the coupler variables

  data LP_NAME / 'LST',       &
                 'SST',       &
                 'SkinT',     &
                 'SkinW',     &
                 'SnowQ',     &
                 'SnowT',     &
                 'XMFLX',     &
                 'YMFLX',     &
                 'ZMFLX',     &
                 'SWUFLX',    &
                 'LWUFLX',    &
                 'SHFLX',     &
                 'LHFLX',     &
                 'QVFLX_Atm', &
                 'GHFLX',     &
                 'PRECFLX',   &
                 'QVFLX_Lnd'  /

  data LP_DESC / 'land surface temp.',               &
                 'sea surface temp.',                &
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
                 'precipitation flux',               &
                 'moisture flux for land'            /

  data LP_UNIT / 'K',       &
                 'K',       &
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
       CPL_LST_UPDATE

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

    allocate( LST  (IA,JA) )
    allocate( SST  (IA,JA) )
    allocate( SkinT(IA,JA) )
    allocate( SkinW(IA,JA) )
    allocate( SnowQ(IA,JA) )
    allocate( SnowT(IA,JA) )

    allocate( XMFLX    (IA,JA) )
    allocate( YMFLX    (IA,JA) )
    allocate( ZMFLX    (IA,JA) )
    allocate( SWUFLX   (IA,JA) )
    allocate( LWUFLX   (IA,JA) )
    allocate( SHFLX    (IA,JA) )
    allocate( LHFLX    (IA,JA) )
    allocate( QVFLX_Atm(IA,JA) )

    allocate( GHFLX    (IA,JA) )
    allocate( PRECFLX  (IA,JA) )
    allocate( QVFLX_Lnd(IA,JA) )

    allocate( DENS(IA,JA) )
    allocate( MOMX(IA,JA) )
    allocate( MOMY(IA,JA) )
    allocate( MOMZ(IA,JA) )
    allocate( RHOS(IA,JA) )
    allocate( PRES(IA,JA) )
    allocate( ATMP(IA,JA) )
    allocate( QV  (IA,JA) )
    allocate( PREC(IA,JA) )
    allocate( SWD (IA,JA) )
    allocate( LWD (IA,JA) )

    allocate( TG  (IA,JA) )
    allocate( QVEF(IA,JA) )
    allocate( EMIT(IA,JA) )
    allocate( ALB (IA,JA) )
    allocate( TCS (IA,JA) )
    allocate( DZG (IA,JA) )
    allocate( Z0M (IA,JA) )
    allocate( Z0H (IA,JA) )
    allocate( Z0E (IA,JA) )

    allocate( Lnd_XMFLX    (IA,JA) )
    allocate( Lnd_YMFLX    (IA,JA) )
    allocate( Lnd_ZMFLX    (IA,JA) )
    allocate( Lnd_SWUFLX   (IA,JA) )
    allocate( Lnd_LWUFLX   (IA,JA) )
    allocate( Lnd_SHFLX    (IA,JA) )
    allocate( Lnd_LHFLX    (IA,JA) )
    allocate( Lnd_QVFLX_Atm(IA,JA) )

    allocate( Lnd_GHFLX    (IA,JA) )
    allocate( Lnd_PRECFLX  (IA,JA) )
    allocate( Lnd_QVFLX_Lnd(IA,JA) )

    Lnd_XMFLX    (:,:) = 0.0_RP
    Lnd_YMFLX    (:,:) = 0.0_RP
    Lnd_ZMFLX    (:,:) = 0.0_RP
    Lnd_SWUFLX   (:,:) = 0.0_RP
    Lnd_LWUFLX   (:,:) = 0.0_RP
    Lnd_SHFLX    (:,:) = 0.0_RP
    Lnd_LHFLX    (:,:) = 0.0_RP
    Lnd_QVFLX_Atm(:,:) = 0.0_RP

    Lnd_GHFLX    (:,:) = 0.0_RP
    Lnd_PRECFLX  (:,:) = 0.0_RP
    Lnd_QVFLX_Lnd(:,:) = 0.0_RP


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
    call COMM_vars8( SkinT(:,:), 3 )
    call COMM_vars8( SkinW(:,:), 4 )
    call COMM_vars8( SnowQ(:,:), 5 )
    call COMM_vars8( SnowT(:,:), 6 )

    call COMM_wait ( LST  (:,:), 1 )
    call COMM_wait ( SST  (:,:), 2 )
    call COMM_wait ( SkinT(:,:), 3 )
    call COMM_wait ( SkinW(:,:), 4 )
    call COMM_wait ( SnowQ(:,:), 5 )
    call COMM_wait ( SnowT(:,:), 6 )

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

!       LST  (:,:) = CONST_UNDEF
       LST  (:,:) = 300.0_RP
       SST  (:,:) = CONST_UNDEF
!       SkinT(:,:) = CONST_UNDEF
       SkinT(:,:) = LST(:,:)
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
       call VALCHECK( SkinT(:,:), 0.0_RP, 1000.0_RP, LP_NAME(I_SkinT), __FILE__, __LINE__ )
       call VALCHECK( SkinW(:,:), 0.0_RP, 1000.0_RP, LP_NAME(I_SkinW), __FILE__, __LINE__ )
       call VALCHECK( SnowQ(:,:), 0.0_RP, 1000.0_RP, LP_NAME(I_SnowQ), __FILE__, __LINE__ )
       call VALCHECK( SnowT(:,:), 0.0_RP, 1000.0_RP, LP_NAME(I_SnowT), __FILE__, __LINE__ )

       call VALCHECK( XMFLX(:,:),     -1.0E4_RP, 1.0E4_RP, LP_NAME(I_XMFLX)    , __FILE__, __LINE__ )
       call VALCHECK( YMFLX(:,:),     -1.0E4_RP, 1.0E4_RP, LP_NAME(I_YMFLX)    , __FILE__, __LINE__ )
       call VALCHECK( ZMFLX(:,:),     -1.0E4_RP, 1.0E4_RP, LP_NAME(I_ZMFLX)    , __FILE__, __LINE__ )
       call VALCHECK( SWUFLX(:,:),    -1.0E4_RP, 1.0E4_RP, LP_NAME(I_SWUFLX)   , __FILE__, __LINE__ )
       call VALCHECK( LWUFLX(:,:),    -1.0E4_RP, 1.0E4_RP, LP_NAME(I_LWUFLX)   , __FILE__, __LINE__ )
       call VALCHECK( SHFLX(:,:),     -1.0E4_RP, 1.0E4_RP, LP_NAME(I_SHFLX)    , __FILE__, __LINE__ )
       call VALCHECK( LHFLX(:,:),     -1.0E4_RP, 1.0E4_RP, LP_NAME(I_LHFLX)    , __FILE__, __LINE__ )
       call VALCHECK( QVFLX_Atm(:,:), -1.0E4_RP, 1.0E4_RP, LP_NAME(I_QVFLX_Atm), __FILE__, __LINE__ )

       call VALCHECK( GHFLX(:,:),     -1.0E4_RP, 1.0E4_RP, LP_NAME(I_GHFLX)    , __FILE__, __LINE__ )
       call VALCHECK( PRECFLX(:,:),   -1.0E4_RP, 1.0E4_RP, LP_NAME(I_PRECFLX)  , __FILE__, __LINE__ )
       call VALCHECK( QVFLX_Lnd(:,:), -1.0E4_RP, 1.0E4_RP, LP_NAME(I_QVFLX_Lnd), __FILE__, __LINE__ )
    endif

    call HIST_in( LST  (:,:), 'LST',   LP_DESC(I_LST),   LP_UNIT(I_LST),   TIME_DTSEC_CPL )
    call HIST_in( SST  (:,:), 'SST',   LP_DESC(I_SST),   LP_UNIT(I_SST),   TIME_DTSEC_CPL )
    call HIST_in( SkinT(:,:), 'SkinT', LP_DESC(I_SkinT), LP_UNIT(I_SkinT), TIME_DTSEC_CPL )
    call HIST_in( SkinW(:,:), 'SkinW', LP_DESC(I_SkinW), LP_UNIT(I_SkinW), TIME_DTSEC_CPL )
    call HIST_in( SnowQ(:,:), 'SnowQ', LP_DESC(I_SnowQ), LP_UNIT(I_SnowQ), TIME_DTSEC_CPL )
    call HIST_in( SnowT(:,:), 'SnowT', LP_DESC(I_SnowT), LP_UNIT(I_SnowT), TIME_DTSEC_CPL )

    call HIST_in( XMFLX    (:,:), 'XMFLX',     LP_DESC(I_XMFLX),     LP_UNIT(I_XMFLX),     TIME_DTSEC_CPL )
    call HIST_in( YMFLX    (:,:), 'YMFLX',     LP_DESC(I_YMFLX),     LP_UNIT(I_YMFLX),     TIME_DTSEC_CPL )
    call HIST_in( ZMFLX    (:,:), 'ZMFLX',     LP_DESC(I_ZMFLX),     LP_UNIT(I_ZMFLX),     TIME_DTSEC_CPL )
    call HIST_in( SWUFLX   (:,:), 'SWUFLX',    LP_DESC(I_SWUFLX),    LP_UNIT(I_SWUFLX),    TIME_DTSEC_CPL )
    call HIST_in( LWUFLX   (:,:), 'LWUFLX',    LP_DESC(I_LWUFLX),    LP_UNIT(I_LWUFLX),    TIME_DTSEC_CPL )
    call HIST_in( SHFLX    (:,:), 'SHFLX',     LP_DESC(I_SHFLX),     LP_UNIT(I_SHFLX),     TIME_DTSEC_CPL )
    call HIST_in( LHFLX    (:,:), 'LHFLX',     LP_DESC(I_LHFLX),     LP_UNIT(I_LHFLX),     TIME_DTSEC_CPL )
    call HIST_in( QVFLX_Atm(:,:), 'QVFLX_Atm', LP_DESC(I_QVFLX_Atm), LP_UNIT(I_QVFLX_Atm), TIME_DTSEC_CPL )

    call HIST_in( GHFLX    (:,:), 'GHFLX',     LP_DESC(I_GHFLX),     LP_UNIT(I_GHFLX),     TIME_DTSEC_CPL )
    call HIST_in( PRECFLX  (:,:), 'PRECFLX',   LP_DESC(I_PRECFLX),   LP_UNIT(I_PRECFLX),   TIME_DTSEC_CPL )
    call HIST_in( QVFLX_Lnd(:,:), 'QVFLX_Lnd', LP_DESC(I_QVFLX_Lnd), LP_UNIT(I_QVFLX_Lnd), TIME_DTSEC_CPL )

    return
  end subroutine CPL_vars_history

  !-----------------------------------------------------------------------------
  !> Budget monitor for coupler
  subroutine CPL_vars_total
!    use mod_comm, only: &
!       STAT_checktotal, &
!       STAT_total
    implicit none

    !real(RP) :: total
    !---------------------------------------------------------------------------

!    if ( STAT_checktotal ) then
!
!       call STAT_total( total, LST(:,:),     LP_NAME(I_LST)   )
!       call STAT_total( total, SST(:,:),     LP_NAME(I_SST)   )
!       call STAT_total( total, SkinT(:,:),   LP_NAME(I_SkinT) )
!       call STAT_total( total, SkinW(:,:),   LP_NAME(I_SkinW) )
!       call STAT_total( total, SnowQ(:,:),   LP_NAME(I_SnowQ) )
!       call STAT_total( total, SnowT(:,:),   LP_NAME(I_SnowT) )
!
!    endif

    return
  end subroutine CPL_vars_total

  subroutine CPL_vars_merge
    implicit none

    ! merge Land-Ocean
    if ( CNT_getCPL2Atm > 0 ) then
       XMFLX    (:,:) = Lnd_XMFLX    (:,:) / CNT_getCPL2Atm
       YMFLX    (:,:) = Lnd_YMFLX    (:,:) / CNT_getCPL2Atm
       ZMFLX    (:,:) = Lnd_ZMFLX    (:,:) / CNT_getCPL2Atm
       SWUFLX   (:,:) = Lnd_SWUFLX   (:,:) / CNT_getCPL2Atm
       LWUFLX   (:,:) = Lnd_LWUFLX   (:,:) / CNT_getCPL2Atm
       SHFLX    (:,:) = Lnd_SHFLX    (:,:) / CNT_getCPL2Atm
       LHFLX    (:,:) = Lnd_LHFLX    (:,:) / CNT_getCPL2Atm
       QVFLX_Atm(:,:) = Lnd_QVFLX_Atm(:,:) / CNT_getCPL2Atm
    end if

    if ( CNT_getCPL2Lnd > 0 ) then
       GHFLX    (:,:) = Lnd_GHFLX    (:,:) / CNT_getCPL2Lnd
       PRECFLX  (:,:) = Lnd_PRECFLX  (:,:) / CNT_getCPL2Lnd
       QVFLX_Lnd(:,:) = Lnd_QVFLX_Lnd(:,:) / CNT_getCPL2Lnd
    end if

    SkinT(:,:) = LST(:,:)

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

    DENS(:,:) = DENS(:,:) + pDENS(:,:)
    MOMX(:,:) = MOMX(:,:) + pMOMX(:,:)
    MOMY(:,:) = MOMY(:,:) + pMOMY(:,:)
    MOMZ(:,:) = MOMZ(:,:) + pMOMZ(:,:)
    RHOS(:,:) = RHOS(:,:) + pRHOS(:,:)
    PRES(:,:) = PRES(:,:) + pPRES(:,:)
    ATMP(:,:) = ATMP(:,:) + pATMP(:,:)
    QV  (:,:) = QV  (:,:) + pQV  (:,:)
    PREC(:,:) = PREC(:,:) + pPREC(:,:)
    SWD (:,:) = SWD (:,:) + pSWD (:,:)
    LWD (:,:) = LWD (:,:) + pLWD (:,:)

    CNT_putAtm = CNT_putAtm + 1.0_RP

    return
  end subroutine CPL_putAtm

  subroutine CPL_putLnd( &
      pTG, pQVEF, pEMIT, & ! (in)
      pALB, pTCS, pDZG,  & ! (in)
      pZ0M, pZ0H, pZ0E   ) ! (in)
    implicit none

    real(RP), intent(in) :: pTG  (IA,JA)
    real(RP), intent(in) :: pQVEF(IA,JA)
    real(RP), intent(in) :: pEMIT(IA,JA)
    real(RP), intent(in) :: pALB (IA,JA)
    real(RP), intent(in) :: pTCS (IA,JA)
    real(RP), intent(in) :: pDZG (IA,JA)
    real(RP), intent(in) :: pZ0M (IA,JA)
    real(RP), intent(in) :: pZ0H (IA,JA)
    real(RP), intent(in) :: pZ0E (IA,JA)

    TG  (:,:) = TG  (:,:) + pTG  (:,:)
    QVEF(:,:) = QVEF(:,:) + pQVEF(:,:)
    EMIT(:,:) = EMIT(:,:) + pEMIT(:,:)
    ALB (:,:) = ALB (:,:) + pALB (:,:)
    TCS (:,:) = TCS (:,:) + pTCS (:,:)
    DZG (:,:) = DZG (:,:) + pDZG (:,:)
    Z0M (:,:) = Z0M (:,:) + pZ0M (:,:)
    Z0H (:,:) = Z0H (:,:) + pZ0H (:,:)
    Z0E (:,:) = Z0E (:,:) + pZ0E (:,:)

    CNT_putLnd = CNT_putLnd + 1.0_RP

    return
  end subroutine CPL_putLnd

  subroutine CPL_AtmLnd_putCPL( &
      pXMFLX, pYMFLX, pZMFLX, & ! (in)
      pSWUFLX, pLWUFLX,       & ! (in)
      pSHFLX, pLHFLX,         & ! (in)
      pGHFLX, pPRECFLX        ) ! (in)
    use mod_const, only: &
       LH0 => CONST_LH0
    implicit none

    real(RP), intent(in) :: pXMFLX  (IA,JA)
    real(RP), intent(in) :: pYMFLX  (IA,JA)
    real(RP), intent(in) :: pZMFLX  (IA,JA)
    real(RP), intent(in) :: pSWUFLX (IA,JA)
    real(RP), intent(in) :: pLWUFLX (IA,JA)
    real(RP), intent(in) :: pSHFLX  (IA,JA)
    real(RP), intent(in) :: pLHFLX  (IA,JA)
    real(RP), intent(in) :: pGHFLX  (IA,JA)
    real(RP), intent(in) :: pPRECFLX(IA,JA)

    Lnd_XMFLX    (:,:) = Lnd_XMFLX    (:,:) + pXMFLX (:,:)
    Lnd_YMFLX    (:,:) = Lnd_YMFLX    (:,:) + pYMFLX (:,:)
    Lnd_ZMFLX    (:,:) = Lnd_ZMFLX    (:,:) + pZMFLX (:,:)
    Lnd_SWUFLX   (:,:) = Lnd_SWUFLX   (:,:) + pSWUFLX(:,:)
    Lnd_LWUFLX   (:,:) = Lnd_LWUFLX   (:,:) + pLWUFLX(:,:)
    Lnd_SHFLX    (:,:) = Lnd_SHFLX    (:,:) + pSHFLX (:,:)
    Lnd_LHFLX    (:,:) = Lnd_LHFLX    (:,:) + pLHFLX (:,:)
    Lnd_QVFLX_Atm(:,:) = Lnd_QVFLX_Atm(:,:) + pLHFLX (:,:)/LH0

    Lnd_GHFLX    (:,:) = Lnd_GHFLX    (:,:) + pGHFLX  (:,:)
    Lnd_PRECFLX  (:,:) = Lnd_PRECFLX  (:,:) + pPRECFLX(:,:)
    Lnd_QVFLX_Lnd(:,:) = Lnd_QVFLX_Lnd(:,:) - pLHFLX  (:,:)/LH0

    CNT_getCPL2Atm = CNT_getCPL2Atm + 1.0_RP
    CNT_getCPL2Lnd = CNT_getCPL2Lnd + 1.0_RP

    return
  end subroutine CPL_AtmLnd_putCPL

  subroutine CPL_getCPL2Atm( &
      pXMFLX, pYMFLX, pZMFLX, & ! (out)
      pSWUFLX, pLWUFLX,       & ! (out)
      pSHFLX, pLHFLX,         & ! (out)
      pQVFLX_Atm              ) ! (out)
    implicit none

    real(RP), intent(out) :: pXMFLX    (IA,JA)
    real(RP), intent(out) :: pYMFLX    (IA,JA)
    real(RP), intent(out) :: pZMFLX    (IA,JA)
    real(RP), intent(out) :: pSWUFLX   (IA,JA)
    real(RP), intent(out) :: pLWUFLX   (IA,JA)
    real(RP), intent(out) :: pSHFLX    (IA,JA)
    real(RP), intent(out) :: pLHFLX    (IA,JA)
    real(RP), intent(out) :: pQVFLX_Atm(IA,JA)

    pXMFLX    (:,:) = XMFLX    (:,:)
    pYMFLX    (:,:) = YMFLX    (:,:)
    pZMFLX    (:,:) = ZMFLX    (:,:)
    pSWUFLX   (:,:) = SWUFLX   (:,:)
    pLWUFLX   (:,:) = LWUFLX   (:,:)
    pSHFLX    (:,:) = SHFLX    (:,:)
    pLHFLX    (:,:) = LHFLX    (:,:)
    pQVFLX_Atm(:,:) = QVFLX_Atm(:,:)

    return
  end subroutine CPL_getCPL2Atm

  subroutine CPL_getCPL2Lnd( &
      pGHFLX, pPRECFLX, pQVFLX_Lnd ) ! (out)
    implicit none

    real(RP), intent(out) :: pGHFLX    (IA,JA)
    real(RP), intent(out) :: pPRECFLX  (IA,JA)
    real(RP), intent(out) :: pQVFLX_Lnd(IA,JA)

    pGHFLX    (:,:) = GHFLX    (:,:)
    pPRECFLX  (:,:) = PRECFLX  (:,:)
    pQVFLX_Lnd(:,:) = QVFLX_Lnd(:,:)

    return
  end subroutine CPL_getCPL2Lnd

  subroutine CPL_AtmLnd_getAtm2CPL( &
      pDENS, pMOMX, pMOMY, pMOMZ, & ! (out)
      pRHOS, pPRES, pATMP, pQV,   & ! (out)
      pPREC, pSWD, pLWD           ) ! (out)
    implicit none

    real(RP), intent(out) :: pDENS(IA,JA)
    real(RP), intent(out) :: pMOMX(IA,JA)
    real(RP), intent(out) :: pMOMY(IA,JA)
    real(RP), intent(out) :: pMOMZ(IA,JA)
    real(RP), intent(out) :: pRHOS(IA,JA)
    real(RP), intent(out) :: pPRES(IA,JA)
    real(RP), intent(out) :: pATMP(IA,JA)
    real(RP), intent(out) :: pQV  (IA,JA)
    real(RP), intent(out) :: pPREC(IA,JA)
    real(RP), intent(out) :: pSWD (IA,JA)
    real(RP), intent(out) :: pLWD (IA,JA)

    pDENS(:,:) = DENS(:,:) / CNT_putAtm
    pMOMX(:,:) = MOMX(:,:) / CNT_putAtm
    pMOMY(:,:) = MOMY(:,:) / CNT_putAtm
    pMOMZ(:,:) = MOMZ(:,:) / CNT_putAtm
    pRHOS(:,:) = RHOS(:,:) / CNT_putAtm
    pPRES(:,:) = PRES(:,:) / CNT_putAtm
    pATMP(:,:) = ATMP(:,:) / CNT_putAtm
    pQV  (:,:) = QV  (:,:) / CNT_putAtm
    pPREC(:,:) = PREC(:,:) / CNT_putAtm
    pSWD (:,:) = SWD (:,:) / CNT_putAtm
    pLWD (:,:) = LWD (:,:) / CNT_putAtm

    return
  end subroutine CPL_AtmLnd_getAtm2CPL

  subroutine CPL_AtmLnd_getLnd2CPL( &
      pTG, pQVEF, pEMIT, & ! (out)
      pALB, pTCS, pDZG,  & ! (out)
      pZ0M, pZ0H, pZ0E   ) ! (out)
    implicit none

    real(RP), intent(out) :: pTG  (IA,JA)
    real(RP), intent(out) :: pQVEF(IA,JA)
    real(RP), intent(out) :: pEMIT(IA,JA)
    real(RP), intent(out) :: pALB (IA,JA)
    real(RP), intent(out) :: pTCS (IA,JA)
    real(RP), intent(out) :: pDZG (IA,JA)
    real(RP), intent(out) :: pZ0M (IA,JA)
    real(RP), intent(out) :: pZ0H (IA,JA)
    real(RP), intent(out) :: pZ0E (IA,JA)

    pTG  (:,:)    = TG  (:,:) / CNT_putLnd
    pQVEF(:,:)    = QVEF(:,:) / CNT_putLnd
    pEMIT(:,:)    = EMIT(:,:) / CNT_putLnd
    pALB (:,:)    = ALB (:,:) / CNT_putLnd
    pTCS (:,:)    = TCS (:,:) / CNT_putLnd
    pDZG (:,:)    = DZG (:,:) / CNT_putLnd
    pZ0M (:,:)    = Z0M (:,:) / CNT_putLnd
    pZ0H (:,:)    = Z0H (:,:) / CNT_putLnd
    pZ0E (:,:)    = Z0E (:,:) / CNT_putLnd

    return
  end subroutine CPL_AtmLnd_getLnd2CPL

  subroutine CPL_flushAtm
    implicit none

    DENS(:,:) = 0.0_RP
    MOMX(:,:) = 0.0_RP
    MOMY(:,:) = 0.0_RP
    MOMZ(:,:) = 0.0_RP
    RHOS(:,:) = 0.0_RP
    PRES(:,:) = 0.0_RP
    ATMP(:,:) = 0.0_RP
    QV  (:,:) = 0.0_RP
    PREC(:,:) = 0.0_RP
    SWD (:,:) = 0.0_RP
    LWD (:,:) = 0.0_RP

    Lnd_XMFLX    (:,:) = 0.0_RP
    Lnd_YMFLX    (:,:) = 0.0_RP
    Lnd_ZMFLX    (:,:) = 0.0_RP
    Lnd_SWUFLX   (:,:) = 0.0_RP
    Lnd_LWUFLX   (:,:) = 0.0_RP
    Lnd_SHFLX    (:,:) = 0.0_RP
    Lnd_LHFLX    (:,:) = 0.0_RP
    Lnd_QVFLX_Atm(:,:) = 0.0_RP

    CNT_getCPL2Atm = 0.0_RP

    return
  end subroutine CPL_flushAtm

  subroutine CPL_flushLnd
    implicit none

    TG  (:,:) = 0.0_RP
    QVEF(:,:) = 0.0_RP
    EMIT(:,:) = 0.0_RP
    ALB (:,:) = 0.0_RP
    TCS (:,:) = 0.0_RP
    DZG (:,:) = 0.0_RP
    Z0M (:,:) = 0.0_RP
    Z0H (:,:) = 0.0_RP
    Z0E (:,:) = 0.0_RP

    Lnd_GHFLX    (:,:) = 0.0_RP
    Lnd_PRECFLX  (:,:) = 0.0_RP
    Lnd_QVFLX_Lnd(:,:) = 0.0_RP

    CNT_getCPL2Lnd = 0.0_RP

    return
  end subroutine CPL_flushLnd

  subroutine CPL_AtmLnd_flushCPL
    implicit none

    CNT_putAtm     = 0.0_RP
    CNT_putLnd     = 0.0_RP

    return
  end subroutine CPL_AtmLnd_flushCPL

end module mod_CPL_vars
