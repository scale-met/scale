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

  character(len=H_SHORT), public, save :: CPL_TYPE_AtmLnd = 'OFF' !< atmos-land coupler type
  character(len=H_SHORT), public, save :: CPL_TYPE_AtmOcn = 'OFF' !< atmos-ocean coupler type
  logical,                public, save :: CPL_sw_AtmLnd           !< do atmos-land coupler calculation?
  logical,                public, save :: CPL_sw_AtmOcn           !< do atmos-ocean coupler calculation?
  logical,                public, save :: CPL_sw_restart          !< output coupler restart?
  character(len=H_SHORT), public, save :: CPL_BULK_TYPE   = ''    !< what is bulk method?

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
  real(RP), private, allocatable :: SFLX_MOMX (:,:) ! momentum flux for x [kg/m2/s]
  real(RP), private, allocatable :: SFLX_MOMY (:,:) ! momentum flux for y [kg/m2/s]
  real(RP), private, allocatable :: SFLX_MOMZ (:,:) ! momentum flux for z [kg/m2/s]
  real(RP), private, allocatable :: SFLX_SWU  (:,:) ! upward short-wave radiation flux (upward positive) [W/m2]
  real(RP), private, allocatable :: SFLX_LWU  (:,:) ! upward long-wave radiation flux (upward positive) [W/m2]
  real(RP), private, allocatable :: SFLX_SH   (:,:) ! sensible heat flux (upward positive) [W/m2]
  real(RP), private, allocatable :: SFLX_LH   (:,:) ! latent heat flux (upward positive) [W/m2]
  real(RP), private, allocatable :: SFLX_QVAtm(:,:) ! moisture flux for atmosphere [kg/m2/s]

  ! surface fluxes for land
  real(RP), private, allocatable :: SFLX_GH   (:,:) ! ground heat flux (upward positive) [W/m2]
  real(RP), private, allocatable :: SFLX_PREC (:,:) ! precipitation flux [kg/m2/s]
  real(RP), private, allocatable :: SFLX_QVLnd(:,:) ! moisture flux for land [kg/m2/s]

  ! Atmospheric values
  real(RP), private, allocatable :: DENS(:,:,:)   ! air density [kg/m3]
  real(RP), private, allocatable :: MOMX(:,:,:)   ! momentum x [kg/m2/s]
  real(RP), private, allocatable :: MOMY(:,:,:)   ! momentum y [kg/m2/s]
  real(RP), private, allocatable :: MOMZ(:,:,:)   ! momentum z [kg/m2/s]
  real(RP), private, allocatable :: RHOT(:,:,:)   ! rho * theta [K*kg/m3]
  real(RP), private, allocatable :: QTRC(:,:,:,:) ! ratio of mass of tracer to total mass [kg/kg]
  real(RP), private, allocatable :: PREC(:,:)       ! surface precipitation rate [kg/m2/s]
  real(RP), private, allocatable :: SWD (:,:)       ! downward short-wave radiation flux (upward positive) [W/m2]
  real(RP), private, allocatable :: LWD (:,:)       ! downward long-wave radiation flux (upward positive) [W/m2]

  ! Land values
  real(RP), private, allocatable :: TG   (:,:) ! soil temperature [K]
  real(RP), private, allocatable :: QvEfc(:,:) ! efficiency of evaporation [no unit]
  real(RP), private, allocatable :: EMIT (:,:) ! emissivity in long-wave radiation [no unit]
  real(RP), private, allocatable :: ALB  (:,:) ! surface albedo in short-wave radiation [no unit]
  real(RP), private, allocatable :: TCS  (:,:) ! thermal conductivity for soil [W/m/K]
  real(RP), private, allocatable :: DZg  (:,:) ! soil depth [m]
  real(RP), private, allocatable :: Z0M  (:,:) ! roughness length for momemtum [m]
  real(RP), private, allocatable :: Z0H  (:,:) ! roughness length for heat [m]
  real(RP), private, allocatable :: Z0E  (:,:) ! roughness length for moisture [m]

  ! AtmLnd surface fluxes for atmosphere
  real(RP), private, allocatable :: Lnd_SFLX_MOMX (:,:) ! momentum flux for x [kg/m2/s]
  real(RP), private, allocatable :: Lnd_SFLX_MOMY (:,:) ! momentum flux for y [kg/m2/s]
  real(RP), private, allocatable :: Lnd_SFLX_MOMZ (:,:) ! momentum flux for z [kg/m2/s]
  real(RP), private, allocatable :: Lnd_SFLX_SWU  (:,:) ! upward short-wave radiation flux (upward positive) [W/m2]
  real(RP), private, allocatable :: Lnd_SFLX_LWU  (:,:) ! upward long-wave radiation flux (upward positive) [W/m2]
  real(RP), private, allocatable :: Lnd_SFLX_SH   (:,:) ! sensible heat flux (upward positive) [W/m2]
  real(RP), private, allocatable :: Lnd_SFLX_LH   (:,:) ! latent heat flux (upward positive) [W/m2]
  real(RP), private, allocatable :: Lnd_SFLX_QVAtm(:,:) ! moisture flux for atmosphere [kg/m2/s]

  ! AtmLnd surface fluxes for land
  real(RP), private, allocatable :: Lnd_SFLX_GH   (:,:) ! ground heat flux (upward positive) [W/m2]
  real(RP), private, allocatable :: Lnd_SFLX_PREC (:,:) ! precipitation flux [kg/m2/s]
  real(RP), private, allocatable :: Lnd_SFLX_QVLnd(:,:) ! moisture flux for land [kg/m2/s]

  ! counter
  real(RP), private, save :: CNT_putAtm     ! counter for putAtm
  real(RP), private, save :: CNT_putLnd     ! counter for putLnd
  real(RP), private, save :: CNT_getCPL2Atm ! counter for getDat2Atm
  real(RP), private, save :: CNT_getCPL2Lnd ! counter for getDat2Lnd

  integer,                    private, save :: I_LST        = 1
  integer,                    private, save :: I_SST        = 2
  integer,                    private, save :: I_SkinT      = 3
  integer,                    private, save :: I_SkinW      = 4
  integer,                    private, save :: I_SnowQ      = 5
  integer,                    private, save :: I_SnowT      = 6
  integer,                    private, save :: I_SFLX_MOMX  = 7
  integer,                    private, save :: I_SFLX_MOMY  = 8
  integer,                    private, save :: I_SFLX_MOMZ  = 9
  integer,                    private, save :: I_SFLX_SWU   = 10
  integer,                    private, save :: I_SFLX_LWU   = 11
  integer,                    private, save :: I_SFLX_SH    = 12
  integer,                    private, save :: I_SFLX_LH    = 13
  integer,                    private, save :: I_SFLX_QVAtm = 14
  integer,                    private, save :: I_SFLX_GH    = 15
  integer,                    private, save :: I_SFLX_PREC  = 16
  integer,                    private, save :: I_SFLX_QVLnd = 17

  character(len=H_SHORT), private, save :: LP_NAME(17) !< name  of the coupler variables
  character(len=H_MID),   private, save :: LP_DESC(17) !< desc. of the coupler variables
  character(len=H_SHORT), private, save :: LP_UNIT(17) !< unit  of the coupler variables

  data LP_NAME / 'LST',        &
                 'SST',        &
                 'SkinT',      &
                 'SkinW',      &
                 'SnowQ',      &
                 'SnowT',      &
                 'SFLX_MOMX',  &
                 'SFLX_MOMY',  &
                 'SFLX_MOMZ',  &
                 'SFLX_SWU',   &
                 'SFLX_LWU',   &
                 'SFLX_SH',    &
                 'SFLX_LH',    &
                 'SFLX_QVAtm', &
                 'SFLX_GH',    &
                 'SFLX_PREC',  &
                 'SFLX_QVLnd'  /

  data LP_DESC / 'land surface temp.',               &
                 'sea surface temp.',                &
                 'ground skin temp.',                &
                 'ground skin water',                &
                 'ground snow amount',               &
                 'ground snow temp.',                &
                 'momentum flux for x',              &
                 'momentum flux for y',              &
                 'momentum flux for z',              &
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
       CPL_TYPE_AtmLnd,   &
       CPL_BULK_TYPE

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

    allocate( SFLX_MOMX (IA,JA) )
    allocate( SFLX_MOMY (IA,JA) )
    allocate( SFLX_MOMZ (IA,JA) )
    allocate( SFLX_SWU  (IA,JA) )
    allocate( SFLX_LWU  (IA,JA) )
    allocate( SFLX_SH   (IA,JA) )
    allocate( SFLX_LH   (IA,JA) )
    allocate( SFLX_QVAtm(IA,JA) )

    allocate( SFLX_GH   (IA,JA) )
    allocate( SFLX_PREC (IA,JA) )
    allocate( SFLX_QVLnd(IA,JA) )

    allocate( DENS(KA,IA,JA)    )
    allocate( MOMX(KA,IA,JA)    )
    allocate( MOMY(KA,IA,JA)    )
    allocate( MOMZ(KA,IA,JA)    )
    allocate( RHOT(KA,IA,JA)    )
    allocate( QTRC(KA,IA,JA,QA) )
    allocate( PREC(IA,JA)       )
    allocate( SWD (IA,JA)       )
    allocate( LWD (IA,JA)       )

    allocate( TG   (IA,JA) )
    allocate( QvEfc(IA,JA) )
    allocate( EMIT (IA,JA) )
    allocate( ALB  (IA,JA) )
    allocate( TCS  (IA,JA) )
    allocate( DZg  (IA,JA) )
    allocate( Z0M  (IA,JA) )
    allocate( Z0H  (IA,JA) )
    allocate( Z0E  (IA,JA) )

    allocate( Lnd_SFLX_MOMX (IA,JA) )
    allocate( Lnd_SFLX_MOMY (IA,JA) )
    allocate( Lnd_SFLX_MOMZ (IA,JA) )
    allocate( Lnd_SFLX_SWU  (IA,JA) )
    allocate( Lnd_SFLX_LWU  (IA,JA) )
    allocate( Lnd_SFLX_SH   (IA,JA) )
    allocate( Lnd_SFLX_LH   (IA,JA) )
    allocate( Lnd_SFLX_QVAtm(IA,JA) )

    allocate( Lnd_SFLX_GH   (IA,JA) )
    allocate( Lnd_SFLX_PREC (IA,JA) )
    allocate( Lnd_SFLX_QVLnd(IA,JA) )

    Lnd_SFLX_MOMX (:,:) = 0.0_RP
    Lnd_SFLX_MOMY (:,:) = 0.0_RP
    Lnd_SFLX_MOMZ (:,:) = 0.0_RP
    Lnd_SFLX_SWU  (:,:) = 0.0_RP
    Lnd_SFLX_LWU  (:,:) = 0.0_RP
    Lnd_SFLX_SH   (:,:) = 0.0_RP
    Lnd_SFLX_LH   (:,:) = 0.0_RP
    Lnd_SFLX_QVAtm(:,:) = 0.0_RP

    Lnd_SFLX_GH   (:,:) = 0.0_RP
    Lnd_SFLX_PREC (:,:) = 0.0_RP
    Lnd_SFLX_QVLnd(:,:) = 0.0_RP


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

       call VALCHECK( SFLX_MOMX(:,:),  -1.0E4_RP, 1.0E4_RP, LP_NAME(I_SFLX_MOMX) , __FILE__, __LINE__ )
       call VALCHECK( SFLX_MOMY(:,:),  -1.0E4_RP, 1.0E4_RP, LP_NAME(I_SFLX_MOMY) , __FILE__, __LINE__ )
       call VALCHECK( SFLX_MOMZ(:,:),  -1.0E4_RP, 1.0E4_RP, LP_NAME(I_SFLX_MOMZ) , __FILE__, __LINE__ )
       call VALCHECK( SFLX_SWU(:,:),   -1.0E4_RP, 1.0E4_RP, LP_NAME(I_SFLX_SWU)  , __FILE__, __LINE__ )
       call VALCHECK( SFLX_LWU(:,:),   -1.0E4_RP, 1.0E4_RP, LP_NAME(I_SFLX_LWU)  , __FILE__, __LINE__ )
       call VALCHECK( SFLX_SH(:,:),    -1.0E4_RP, 1.0E4_RP, LP_NAME(I_SFLX_SH)   , __FILE__, __LINE__ )
       call VALCHECK( SFLX_LH(:,:),    -1.0E4_RP, 1.0E4_RP, LP_NAME(I_SFLX_LH)   , __FILE__, __LINE__ )
       call VALCHECK( SFLX_QVAtm(:,:), -1.0E4_RP, 1.0E4_RP, LP_NAME(I_SFLX_QVAtm), __FILE__, __LINE__ )

       call VALCHECK( SFLX_GH(:,:),    -1.0E4_RP, 1.0E4_RP, LP_NAME(I_SFLX_GH)   , __FILE__, __LINE__ )
       call VALCHECK( SFLX_PREC(:,:),  -1.0E4_RP, 1.0E4_RP, LP_NAME(I_SFLX_PREC) , __FILE__, __LINE__ )
       call VALCHECK( SFLX_QVLnd(:,:), -1.0E4_RP, 1.0E4_RP, LP_NAME(I_SFLX_QVLnd), __FILE__, __LINE__ )
    endif

    call HIST_in( LST  (:,:), 'LST',   LP_DESC(I_LST),   LP_UNIT(I_LST),   TIME_DTSEC_CPL )
    call HIST_in( SST  (:,:), 'SST',   LP_DESC(I_SST),   LP_UNIT(I_SST),   TIME_DTSEC_CPL )
    call HIST_in( SkinT(:,:), 'SkinT', LP_DESC(I_SkinT), LP_UNIT(I_SkinT), TIME_DTSEC_CPL )
    call HIST_in( SkinW(:,:), 'SkinW', LP_DESC(I_SkinW), LP_UNIT(I_SkinW), TIME_DTSEC_CPL )
    call HIST_in( SnowQ(:,:), 'SnowQ', LP_DESC(I_SnowQ), LP_UNIT(I_SnowQ), TIME_DTSEC_CPL )
    call HIST_in( SnowT(:,:), 'SnowT', LP_DESC(I_SnowT), LP_UNIT(I_SnowT), TIME_DTSEC_CPL )

    call HIST_in( SFLX_MOMX (:,:), 'SFLX_MOMX',  LP_DESC(I_SFLX_MOMX),  LP_UNIT(I_SFLX_MOMX),  TIME_DTSEC_CPL )
    call HIST_in( SFLX_MOMY (:,:), 'SFLX_MOMY',  LP_DESC(I_SFLX_MOMY),  LP_UNIT(I_SFLX_MOMY),  TIME_DTSEC_CPL )
    call HIST_in( SFLX_MOMZ (:,:), 'SFLX_MOMZ',  LP_DESC(I_SFLX_MOMZ),  LP_UNIT(I_SFLX_MOMZ),  TIME_DTSEC_CPL )
    call HIST_in( SFLX_SWU  (:,:), 'SFLX_SWU',   LP_DESC(I_SFLX_SWU),   LP_UNIT(I_SFLX_SWU),   TIME_DTSEC_CPL )
    call HIST_in( SFLX_LWU  (:,:), 'SFLX_LWU',   LP_DESC(I_SFLX_LWU),   LP_UNIT(I_SFLX_LWU),   TIME_DTSEC_CPL )
    call HIST_in( SFLX_SH   (:,:), 'SFLX_SH',    LP_DESC(I_SFLX_SH),    LP_UNIT(I_SFLX_SH),    TIME_DTSEC_CPL )
    call HIST_in( SFLX_LH   (:,:), 'SFLX_LH',    LP_DESC(I_SFLX_LH),    LP_UNIT(I_SFLX_LH),    TIME_DTSEC_CPL )
    call HIST_in( SFLX_QVAtm(:,:), 'SFLX_QVAtm', LP_DESC(I_SFLX_QVAtm), LP_UNIT(I_SFLX_QVAtm), TIME_DTSEC_CPL )

    call HIST_in( SFLX_GH   (:,:), 'SFLX_GH',    LP_DESC(I_SFLX_GH),    LP_UNIT(I_SFLX_GH),    TIME_DTSEC_CPL )
    call HIST_in( SFLX_PREC (:,:), 'SFLX_PREC',  LP_DESC(I_SFLX_PREC),  LP_UNIT(I_SFLX_PREC),  TIME_DTSEC_CPL )
    call HIST_in( SFLX_QVLnd(:,:), 'SFLX_QVLnd', LP_DESC(I_SFLX_QVLnd), LP_UNIT(I_SFLX_QVLnd), TIME_DTSEC_CPL )

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
       SFLX_MOMX (:,:) = Lnd_SFLX_MOMX (:,:) / CNT_getCPL2Atm
       SFLX_MOMY (:,:) = Lnd_SFLX_MOMY (:,:) / CNT_getCPL2Atm
       SFLX_MOMZ (:,:) = Lnd_SFLX_MOMZ (:,:) / CNT_getCPL2Atm
       SFLX_SWU  (:,:) = Lnd_SFLX_SWU  (:,:) / CNT_getCPL2Atm
       SFLX_LWU  (:,:) = Lnd_SFLX_LWU  (:,:) / CNT_getCPL2Atm
       SFLX_SH   (:,:) = Lnd_SFLX_SH   (:,:) / CNT_getCPL2Atm
       SFLX_LH   (:,:) = Lnd_SFLX_LH   (:,:) / CNT_getCPL2Atm
       SFLX_QVAtm(:,:) = Lnd_SFLX_QVAtm(:,:) / CNT_getCPL2Atm
    end if

    if ( CNT_getCPL2Lnd > 0 ) then
       SFLX_GH   (:,:) = Lnd_SFLX_GH   (:,:) / CNT_getCPL2Lnd
       SFLX_PREC (:,:) = Lnd_SFLX_PREC (:,:) / CNT_getCPL2Lnd
       SFLX_QVLnd(:,:) = Lnd_SFLX_QVLnd(:,:) / CNT_getCPL2Lnd
    end if

    SkinT(:,:) = LST(:,:)

  end subroutine CPL_vars_merge

  subroutine CPL_putAtm( &
      pDENS, pMOMX, pMOMY, pMOMZ, pRHOT, & ! (in)
      pQTRC, pPREC, pSWD, pLWD           ) ! (in)
    implicit none

    real(RP), intent(in) :: pDENS(KA,IA,JA)
    real(RP), intent(in) :: pMOMX(KA,IA,JA)
    real(RP), intent(in) :: pMOMY(KA,IA,JA)
    real(RP), intent(in) :: pMOMZ(KA,IA,JA)
    real(RP), intent(in) :: pRHOT(KA,IA,JA)
    real(RP), intent(in) :: pQTRC(KA,IA,JA,QA)
    real(RP), intent(in) :: pPREC(IA,JA)
    real(RP), intent(in) :: pSWD (IA,JA)
    real(RP), intent(in) :: pLWD (IA,JA)

    DENS(:,:,:)   = DENS(:,:,:)   + pDENS  (:,:,:)
    MOMX(:,:,:)   = MOMX(:,:,:)   + pMOMX  (:,:,:)
    MOMY(:,:,:)   = MOMY(:,:,:)   + pMOMY  (:,:,:)
    MOMZ(:,:,:)   = MOMZ(:,:,:)   + pMOMZ  (:,:,:)
    RHOT(:,:,:)   = RHOT(:,:,:)   + pRHOT  (:,:,:)
    QTRC(:,:,:,:) = QTRC(:,:,:,:) + pQTRC  (:,:,:,:)
    PREC(:,:)     = PREC(:,:)     + pPREC  (:,:)
    SWD (:,:)     = SWD (:,:)     + pSWD(:,:)
    LWD (:,:)     = LWD (:,:)     + pLWD(:,:)

    CNT_putAtm = CNT_putAtm + 1.0_RP

    return
  end subroutine CPL_putAtm

  subroutine CPL_putLnd( &
      pTG, pQvEfc, pEMIT, & ! (in)
      pALB, pTCS, pDZg,   & ! (in)
      pZ0M, pZ0H, pZ0E    ) ! (in)
    implicit none

    real(RP), intent(in) :: pTG   (IA,JA)
    real(RP), intent(in) :: pQvEfc(IA,JA)
    real(RP), intent(in) :: pEMIT (IA,JA)
    real(RP), intent(in) :: pALB  (IA,JA)
    real(RP), intent(in) :: pTCS  (IA,JA)
    real(RP), intent(in) :: pDZg  (IA,JA)
    real(RP), intent(in) :: pZ0M  (IA,JA)
    real(RP), intent(in) :: pZ0H  (IA,JA)
    real(RP), intent(in) :: pZ0E  (IA,JA)

    TG   (:,:) = TG   (:,:) + pTG   (:,:)
    QvEfc(:,:) = QvEfc(:,:) + pQvEfc(:,:)
    EMIT (:,:) = EMIT (:,:) + pEMIT (:,:)
    ALB  (:,:) = ALB  (:,:) + pALB  (:,:)
    TCS  (:,:) = TCS  (:,:) + pTCS  (:,:)
    DZg  (:,:) = DZg  (:,:) + pDZg  (:,:)
    Z0M  (:,:) = Z0M  (:,:) + pZ0M  (:,:)
    Z0H  (:,:) = Z0H  (:,:) + pZ0H  (:,:)
    Z0E  (:,:) = Z0E  (:,:) + pZ0E  (:,:)

    CNT_putLnd = CNT_putLnd + 1.0_RP

    return
  end subroutine CPL_putLnd

  subroutine CPL_AtmLnd_putCPL( &
      pSFLX_MOMX, pSFLX_MOMY, pSFLX_MOMZ,       & ! (in)
      pSFLX_SWU, pSFLX_LWU, pSFLX_SH, pSFLX_LH, & ! (in)
      pSFLX_GH, pSFLX_PREC                      ) ! (in)
    use mod_const, only: &
       LH0 => CONST_LH0
    implicit none

    real(RP), intent(in) :: pSFLX_MOMX(IA,JA)
    real(RP), intent(in) :: pSFLX_MOMY(IA,JA)
    real(RP), intent(in) :: pSFLX_MOMZ(IA,JA)
    real(RP), intent(in) :: pSFLX_SWU (IA,JA)
    real(RP), intent(in) :: pSFLX_LWU (IA,JA)
    real(RP), intent(in) :: pSFLX_SH  (IA,JA)
    real(RP), intent(in) :: pSFLX_LH  (IA,JA)
    real(RP), intent(in) :: pSFLX_GH  (IA,JA)
    real(RP), intent(in) :: pSFLX_PREC(IA,JA)

    Lnd_SFLX_MOMX (:,:) = Lnd_SFLX_MOMX (:,:) + pSFLX_MOMX(:,:)
    Lnd_SFLX_MOMY (:,:) = Lnd_SFLX_MOMY (:,:) + pSFLX_MOMY(:,:)
    Lnd_SFLX_MOMZ (:,:) = Lnd_SFLX_MOMZ (:,:) + pSFLX_MOMZ(:,:)
    Lnd_SFLX_SWU  (:,:) = Lnd_SFLX_SWU  (:,:) + pSFLX_SWU (:,:)
    Lnd_SFLX_LWU  (:,:) = Lnd_SFLX_LWU  (:,:) + pSFLX_LWU (:,:)
    Lnd_SFLX_SH   (:,:) = Lnd_SFLX_SH   (:,:) + pSFLX_SH  (:,:)
    Lnd_SFLX_LH   (:,:) = Lnd_SFLX_LH   (:,:) + pSFLX_LH  (:,:)
    Lnd_SFLX_QVAtm(:,:) = Lnd_SFLX_QVAtm(:,:) + pSFLX_LH  (:,:)/LH0

    Lnd_SFLX_GH   (:,:) = Lnd_SFLX_GH   (:,:) + pSFLX_GH  (:,:)
    Lnd_SFLX_PREC (:,:) = Lnd_SFLX_PREC (:,:) + pSFLX_PREC(:,:)
    Lnd_SFLX_QVLnd(:,:) = Lnd_SFLX_QVLnd(:,:) - pSFLX_LH  (:,:)/LH0

    CNT_getCPL2Atm = CNT_getCPL2Atm + 1.0_RP
    CNT_getCPL2Lnd = CNT_getCPL2Lnd + 1.0_RP

    return
  end subroutine CPL_AtmLnd_putCPL

  subroutine CPL_getCPL2Atm( &
      pSFLX_MOMX, pSFLX_MOMY, pSFLX_MOMZ, pSFLX_SWU, pSFLX_LWU, & ! (out)
      pSFLX_SH, pSFLX_LH, pSFLX_QVAtm                           ) ! (out)
    implicit none

    real(RP), intent(out) :: pSFLX_MOMX (IA,JA)
    real(RP), intent(out) :: pSFLX_MOMY (IA,JA)
    real(RP), intent(out) :: pSFLX_MOMZ (IA,JA)
    real(RP), intent(out) :: pSFLX_SWU  (IA,JA)
    real(RP), intent(out) :: pSFLX_LWU  (IA,JA)
    real(RP), intent(out) :: pSFLX_SH   (IA,JA)
    real(RP), intent(out) :: pSFLX_LH   (IA,JA)
    real(RP), intent(out) :: pSFLX_QVAtm(IA,JA)

    pSFLX_MOMX (:,:) = SFLX_MOMX (:,:)
    pSFLX_MOMY (:,:) = SFLX_MOMY (:,:)
    pSFLX_MOMZ (:,:) = SFLX_MOMZ (:,:)
    pSFLX_SWU  (:,:) = SFLX_SWU  (:,:)
    pSFLX_LWU  (:,:) = SFLX_LWU  (:,:)
    pSFLX_SH   (:,:) = SFLX_SH   (:,:)
    pSFLX_LH   (:,:) = SFLX_LH   (:,:)
    pSFLX_QVAtm(:,:) = SFLX_QVAtm(:,:)

    return
  end subroutine CPL_getCPL2Atm

  subroutine CPL_getCPL2Lnd( &
      pSFLX_GH, pSFLX_PREC, pSFLX_QVLnd ) ! (out)
    implicit none

    real(RP), intent(out) :: pSFLX_GH   (IA,JA)
    real(RP), intent(out) :: pSFLX_PREC (IA,JA)
    real(RP), intent(out) :: pSFLX_QVLnd(IA,JA)

    pSFLX_GH   (:,:) = SFLX_GH   (:,:)
    pSFLX_PREC (:,:) = SFLX_PREC (:,:)
    pSFLX_QVLnd(:,:) = SFLX_QVLnd(:,:)

    return
  end subroutine CPL_getCPL2Lnd

  subroutine CPL_AtmLnd_getAtm2CPL( &
      pDENS, pMOMX, pMOMY, pMOMZ, pRHOT, & ! (out)
      pQTRC, pPREC, pSWD, pLWD           ) ! (out)
    implicit none

    real(RP), intent(out) :: pDENS(KA,IA,JA)
    real(RP), intent(out) :: pMOMX(KA,IA,JA)
    real(RP), intent(out) :: pMOMY(KA,IA,JA)
    real(RP), intent(out) :: pMOMZ(KA,IA,JA)
    real(RP), intent(out) :: pRHOT(KA,IA,JA)
    real(RP), intent(out) :: pQTRC(KA,IA,JA,QA)
    real(RP), intent(out) :: pPREC(IA,JA)
    real(RP), intent(out) :: pSWD (IA,JA)
    real(RP), intent(out) :: pLWD (IA,JA)

    pDENS(:,:,:)   = DENS(:,:,:)   / CNT_putAtm
    pMOMX(:,:,:)   = MOMX(:,:,:)   / CNT_putAtm
    pMOMY(:,:,:)   = MOMY(:,:,:)   / CNT_putAtm
    pMOMZ(:,:,:)   = MOMZ(:,:,:)   / CNT_putAtm
    pRHOT(:,:,:)   = RHOT(:,:,:)   / CNT_putAtm
    pQTRC(:,:,:,:) = QTRC(:,:,:,:) / CNT_putAtm
    pPREC(:,:)     = PREC(:,:)     / CNT_putAtm
    pSWD(:,:)      = SWD(:,:)      / CNT_putAtm
    pLWD(:,:)      = LWD(:,:)      / CNT_putAtm

    return
  end subroutine CPL_AtmLnd_getAtm2CPL

  subroutine CPL_AtmLnd_getLnd2CPL( &
      pTG, pQvEfc, pEMIT, & ! (out)
      pALB, pTCS, pDZg,   & ! (out)
      pZ0M, pZ0H, pZ0E    ) ! (out)
    implicit none

    real(RP), intent(out) :: pTG   (IA,JA)
    real(RP), intent(out) :: pQvEfc(IA,JA)
    real(RP), intent(out) :: pEMIT (IA,JA)
    real(RP), intent(out) :: pALB  (IA,JA)
    real(RP), intent(out) :: pTCS  (IA,JA)
    real(RP), intent(out) :: pDZg  (IA,JA)
    real(RP), intent(out) :: pZ0M  (IA,JA)
    real(RP), intent(out) :: pZ0H  (IA,JA)
    real(RP), intent(out) :: pZ0E  (IA,JA)

    pTG   (:,:)    = TG   (:,:) / CNT_putLnd
    pQvEfc(:,:)    = QvEfc(:,:) / CNT_putLnd
    pEMIT (:,:)    = EMIT (:,:) / CNT_putLnd
    pALB  (:,:)    = ALB  (:,:) / CNT_putLnd
    pTCS  (:,:)    = TCS  (:,:) / CNT_putLnd
    pDZg  (:,:)    = DZg  (:,:) / CNT_putLnd
    pZ0M  (:,:)    = Z0M  (:,:) / CNT_putLnd
    pZ0H  (:,:)    = Z0H  (:,:) / CNT_putLnd
    pZ0E  (:,:)    = Z0E  (:,:) / CNT_putLnd

    return
  end subroutine CPL_AtmLnd_getLnd2CPL

  subroutine CPL_flushAtm
    implicit none

    DENS(:,:,:)    = 0.0_RP
    MOMX(:,:,:)    = 0.0_RP
    MOMY(:,:,:)    = 0.0_RP
    MOMZ(:,:,:)    = 0.0_RP
    RHOT(:,:,:)    = 0.0_RP
    QTRC(:,:,:,:)  = 0.0_RP
    PREC(:,:)      = 0.0_RP
    SWD (:,:)      = 0.0_RP
    LWD (:,:)      = 0.0_RP

    Lnd_SFLX_MOMX (:,:) = 0.0_RP
    Lnd_SFLX_MOMY (:,:) = 0.0_RP
    Lnd_SFLX_MOMZ (:,:) = 0.0_RP
    Lnd_SFLX_SWU  (:,:) = 0.0_RP
    Lnd_SFLX_LWU  (:,:) = 0.0_RP
    Lnd_SFLX_SH   (:,:) = 0.0_RP
    Lnd_SFLX_LH   (:,:) = 0.0_RP
    Lnd_SFLX_QVAtm(:,:) = 0.0_RP

    CNT_getCPL2Atm = 0.0_RP

    return
  end subroutine CPL_flushAtm

  subroutine CPL_flushLnd
    implicit none

    TG   (:,:) = 0.0_RP
    QvEfc(:,:) = 0.0_RP
    EMIT (:,:) = 0.0_RP
    ALB  (:,:) = 0.0_RP
    TCS  (:,:) = 0.0_RP
    DZg  (:,:) = 0.0_RP
    Z0M  (:,:) = 0.0_RP
    Z0H  (:,:) = 0.0_RP
    Z0E  (:,:) = 0.0_RP

    Lnd_SFLX_GH   (:,:) = 0.0_RP
    Lnd_SFLX_PREC (:,:) = 0.0_RP
    Lnd_SFLX_QVLnd(:,:) = 0.0_RP

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
