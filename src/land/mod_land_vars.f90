!-------------------------------------------------------------------------------
!> module LAND Variables
!!
!! @par Description
!!          Container for land variables
!!
!! @author Team SCALE
!! @li      2013-08-31 (T.Yamaura)  [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_land_vars
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use gtool_file_h, only: &
     File_HSHORT, &
     File_HMID,   &
     File_HLONG
  use mod_stdio, only: &
     IO_FID_LOG, &
     IO_L,       &
     IO_SYSCHR,  &
     IO_FILECHR
  use mod_time, only: &
     TIME_rapstart, &
     TIME_rapend
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include "inc_precision.h"
  include 'inc_index.h'
!  include 'inc_land.h'

  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: LAND_vars_setup
  public :: LAND_vars_fillhalo
  public :: LAND_vars_restart_read
  public :: LAND_vars_restart_write
  public :: LAND_vars_history
  public :: LAND_vars_total

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, save :: SFLX_GH  ( IA,JA) ! ground heat flux (upward positive) [W/m2]
  real(RP), public, save :: SFLX_PREC (IA,JA) ! precipitation flux [kg/m2/s]
  real(RP), public, save :: SFLX_QVLnd(IA,JA) ! moisture flux [kg/m2/s]

  real(RP), public, save :: TG   (IA,JA) ! soil temperature [K]
  real(RP), public, save :: QvEfc(IA,JA) ! efficiency of evaporation [no unit]
  real(RP), public, save :: EMIT (IA,JA) ! emissivity in long-wave radiation [no unit]
  real(RP), public, save :: ALB  (IA,JA) ! surface albedo in short-wave radiation [no unit]
  real(RP), public, save :: TCS  (IA,JA) ! thermal conductivity for soil [W/m/K]
  real(RP), public, save :: HCS  (IA,JA) ! heat capacity for soil [J/K]
  real(RP), public, save :: DZg  (IA,JA) ! soil depth [m]
  real(RP), public, save :: Z00  (IA,JA) ! basic factor for momemtum
  real(RP), public, save :: Z0R  (IA,JA) ! rough factor for momemtum
  real(RP), public, save :: Z0S  (IA,JA) ! smooth factor for momemtum
  real(RP), public, save :: Zt0  (IA,JA) ! basic factor for heat
  real(RP), public, save :: ZtR  (IA,JA) ! rough factor for heat
  real(RP), public, save :: ZtS  (IA,JA) ! smooth factor for heat
  real(RP), public, save :: Ze0  (IA,JA) ! basic factor for moisture
  real(RP), public, save :: ZeR  (IA,JA) ! rough factor for moisture
  real(RP), public, save :: ZeS  (IA,JA) ! smooth factor for moisture

!  real(RP), public, save :: SoilT(KA_soil,IA,JA) ! Soil temperature             [K]
!  real(RP), public, save :: SoilW(KA_soil,IA,JA) ! Soil moisture (liquid water) [m3/m3]
!  real(RP), public, save :: SoilI(KA_soil,IA,JA) ! Soil ice      (ice    water) [m3/m3]

  integer,                    public, save :: I_TG    = 1
  integer,                    public, save :: I_QvEfc = 2
  integer,                    public, save :: I_EMIT  = 3
  integer,                    public, save :: I_ALB   = 4
  integer,                    public, save :: I_TCS   = 5
  integer,                    public, save :: I_HCS   = 6
  integer,                    public, save :: I_DZg   = 7
  integer,                    public, save :: I_Z00   = 8
  integer,                    public, save :: I_Z0R   = 9
  integer,                    public, save :: I_Z0S   = 10
  integer,                    public, save :: I_Zt0   = 11
  integer,                    public, save :: I_ZtR   = 12
  integer,                    public, save :: I_ZtS   = 13
  integer,                    public, save :: I_Ze0   = 14
  integer,                    public, save :: I_ZeR   = 15
  integer,                    public, save :: I_ZeS   = 16
  integer,                    public, save :: I_SoilT = 17
  integer,                    public, save :: I_SoilW = 18
  integer,                    public, save :: I_SoilI = 19
  character(len=File_HSHORT), public, save :: LP_NAME(19) !< name  of the land variables
  character(len=File_HMID),   public, save :: LP_DESC(19) !< desc. of the land variables
  character(len=File_HSHORT), public, save :: LP_UNIT(19) !< unit  of the land variables

  data LP_NAME / 'TG',    &
                 'QvEfc', &
                 'EMIT',  &
                 'ALB',   &
                 'TCS',   &
                 'HCS',   &
                 'DZg',   &
                 'Z00',   &
                 'Z0R',   &
                 'Z0S',   &
                 'Zt0',   &
                 'ZtR',   &
                 'ZtS',   &
                 'Ze0',   &
                 'ZeR',   &
                 'ZeS',   &
                 'SoilT', &
                 'SoilW', &
                 'SoilI'  /
  data LP_DESC / 'soil temperature',                       &
                 'efficiency of evaporation',              &
                 'emissivity in long-wave radiation',      &
                 'surface albedo in short-wave radiation', &
                 'thermal conductivity for soil',          &
                 'heat capacity for soil',                 &
                 'soil depth',                             &
                 'basic factor for momemtum',              &
                 'rough factor for momemtum',              &
                 'smooth factor for momemtum',             &
                 'basic factor for heat',                  &
                 'rough factor for heat',                  &
                 'smooth factor for heat',                 &
                 'basic factor for moisture',              &
                 'rough factor for moisture',              &
                 'smooth factor for moisture',             &
                 'soil temperature',                       &
                 'soil moisture',                          &
                 'soil ice'                                /
  data LP_UNIT / 'K',       &
                 'no-unit', &
                 'no-unit', &
                 'no-unit', &
                 'W/m/K',   &
                 'J/K',     &
                 'm',       &
                 'no-unit', &
                 'no-unit', &
                 'no-unit', &
                 'no-unit', &
                 'no-unit', &
                 'no-unit', &
                 'no-unit', &
                 'no-unit', &
                 'no-unit', &
                 'K',       &
                 'm3/m3',   &
                 'm3/m3'    /

  character(len=IO_SYSCHR),  public, save :: LAND_TYPE_PHY = 'OFF' !< Land physics type
  logical,                   public, save :: LAND_sw_phy           !< do land physics update?
  logical,                   public, save :: LAND_sw_restart       !< output restart?

  character(len=IO_FILECHR), public, save :: LAND_RESTART_IN_BASENAME = '' !< basename of the input file

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical,                   private, save :: LAND_RESTART_OUTPUT       = .false.             !< output restart file?
  character(len=IO_FILECHR), private, save :: LAND_RESTART_OUT_BASENAME = 'restart_out'       !< basename of the output file
  character(len=IO_SYSCHR),  private, save :: LAND_RESTART_OUT_TITLE    = 'SCALE3 LAND VARS.' !< title    of the output file
  character(len=IO_SYSCHR),  private, save :: LAND_RESTART_OUT_DTYPE    = 'DEFAULT'           !< REAL4 or REAL8

  logical,                   private, save :: LAND_VARS_CHECKRANGE      = .false.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine LAND_vars_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    NAMELIST / PARAM_LAND / &
       LAND_TYPE_PHY

    NAMELIST / PARAM_LAND_VARS /  &
       LAND_RESTART_IN_BASENAME,  &
       LAND_RESTART_OUTPUT,       &
       LAND_RESTART_OUT_BASENAME, &
       LAND_RESTART_OUT_TITLE,    &
       LAND_VARS_CHECKRANGE

    integer :: ierr
    integer :: ip
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[LAND VARS]/Categ[LAND]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_LAND,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_LAND. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_LAND)

    if( IO_L ) write(IO_FID_LOG,*) '*** [LAND] selected components'

    if ( LAND_TYPE_PHY /= 'OFF' .AND. LAND_TYPE_PHY /= 'NONE' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Land physics : ON'
       LAND_sw_phy = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Land physics : OFF'
       LAND_sw_phy = .false.
    endif

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[LAND VARS]/Categ[LAND]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_LAND_VARS,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_LAND_VARS. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_LAND_VARS)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** [LAND ] prognostic variables'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A8,A,A32,3(A))') &
               '***       |',' VARNAME','|', 'DESCRIPTION                     ','[', 'UNIT            ',']'
    do ip = 1, 19
       if( IO_L ) write(IO_FID_LOG,'(1x,A,i3,A,A8,A,A32,3(A))') &
                  '*** NO.',ip,'|',trim(LP_NAME(ip)),'|', LP_DESC(ip),'[', LP_UNIT(ip),']'
    enddo

    if( IO_L ) write(IO_FID_LOG,*) 'Output...'
    if ( LAND_RESTART_OUTPUT ) then
       if( IO_L ) write(IO_FID_LOG,*) '  Land restart output : YES'
       LAND_sw_restart = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '  Land restart output : NO'
       LAND_sw_restart = .false.
    endif
    if( IO_L ) write(IO_FID_LOG,*)

    return
  end subroutine LAND_vars_setup

  !-----------------------------------------------------------------------------
  !> fill HALO region of land variables
  subroutine LAND_vars_fillhalo
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none
    !---------------------------------------------------------------------------

    ! fill IHALO & JHALO
    call COMM_vars8( TG   (:,:),   1  )
    call COMM_vars8( QvEfc(:,:),   2  )
    call COMM_vars8( EMIT (:,:),   3  )
    call COMM_vars8( ALB  (:,:),   4  )
    call COMM_vars8( TCS  (:,:),   5  )
    call COMM_vars8( HCS  (:,:),   6  )
    call COMM_vars8( DZg  (:,:),   7  )
    call COMM_vars8( Z00  (:,:),   8  )
    call COMM_vars8( Z0R  (:,:),   9  )
    call COMM_vars8( Z0S  (:,:),   10 )
    call COMM_vars8( Zt0  (:,:),   11 )
    call COMM_vars8( ZtR  (:,:),   12 )
    call COMM_vars8( ZtS  (:,:),   13 )
    call COMM_vars8( Ze0  (:,:),   14 )
    call COMM_vars8( ZeR  (:,:),   15 )
    call COMM_vars8( ZeS  (:,:),   16 )

    call COMM_wait ( TG   (:,:),   1  )
    call COMM_wait ( QvEfc(:,:),   2  )
    call COMM_wait ( EMIT (:,:),   3  )
    call COMM_wait ( ALB  (:,:),   4  )
    call COMM_wait ( TCS  (:,:),   5  )
    call COMM_wait ( HCS  (:,:),   6  )
    call COMM_wait ( DZg  (:,:),   7  )
    call COMM_wait ( Z00  (:,:),   8  )
    call COMM_wait ( Z0R  (:,:),   9  )
    call COMM_wait ( Z0S  (:,:),   10 )
    call COMM_wait ( Zt0  (:,:),   11 )
    call COMM_wait ( ZtR  (:,:),   12 )
    call COMM_wait ( ZtS  (:,:),   13 )
    call COMM_wait ( Ze0  (:,:),   14 )
    call COMM_wait ( ZeR  (:,:),   15 )
    call COMM_wait ( ZeS  (:,:),   16 )

!    call COMM_vars8( SoilT(:,:,:), 1 )
!    call COMM_vars8( SoilW(:,:,:), 2 )
!    call COMM_vars8( SoilI(:,:,:), 3 )

!    call COMM_wait ( SoilT(:,:,:), 1 )
!    call COMM_wait ( SoilW(:,:,:), 2 )
!    call COMM_wait ( SoilI(:,:,:), 3 )

    return
  end subroutine LAND_vars_fillhalo

  !-----------------------------------------------------------------------------
  !> Read land restart
  subroutine LAND_vars_restart_read
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
    if( IO_L ) write(IO_FID_LOG,*) '*** Input restart file (land) ***'

    call TIME_rapstart('FILE I NetCDF')

    if ( LAND_RESTART_IN_BASENAME /= '' ) then

       call FILEIO_read( TG(:,:),                                         & ! [OUT]
                         LAND_RESTART_IN_BASENAME, 'TG', 'XY', step=1     ) ! [IN]
       call FILEIO_read( QvEfc(:,:),                                      & ! [OUT]
                         LAND_RESTART_IN_BASENAME, 'QvEfc', 'XY', step=1  ) ! [IN]
       call FILEIO_read( EMIT(:,:),                                       & ! [OUT]
                         LAND_RESTART_IN_BASENAME, 'EMIT', 'XY', step=1   ) ! [IN]
       call FILEIO_read( ALB(:,:),                                        & ! [OUT]
                         LAND_RESTART_IN_BASENAME, 'ALB', 'XY', step=1    ) ! [IN]
       call FILEIO_read( TCS(:,:),                                        & ! [OUT]
                         LAND_RESTART_IN_BASENAME, 'TCS', 'XY', step=1    ) ! [IN]
       call FILEIO_read( HCS(:,:),                                        & ! [OUT]
                         LAND_RESTART_IN_BASENAME, 'HCS', 'XY', step=1    ) ! [IN]
       call FILEIO_read( DZg(:,:),                                        & ! [OUT]
                         LAND_RESTART_IN_BASENAME, 'DZg', 'XY', step=1    ) ! [IN]
       call FILEIO_read( Z00(:,:),                                        & ! [OUT]
                         LAND_RESTART_IN_BASENAME, 'Z00', 'XY', step=1    ) ! [IN]
       call FILEIO_read( Z0R(:,:),                                        & ! [OUT]
                         LAND_RESTART_IN_BASENAME, 'Z0R', 'XY', step=1    ) ! [IN]
       call FILEIO_read( Z0S(:,:),                                        & ! [OUT]
                         LAND_RESTART_IN_BASENAME, 'Z0S', 'XY', step=1    ) ! [IN]
       call FILEIO_read( Zt0(:,:),                                        & ! [OUT]
                         LAND_RESTART_IN_BASENAME, 'Zt0', 'XY', step=1    ) ! [IN]
       call FILEIO_read( ZtR(:,:),                                        & ! [OUT]
                         LAND_RESTART_IN_BASENAME, 'ZtR', 'XY', step=1    ) ! [IN]
       call FILEIO_read( ZtS(:,:),                                        & ! [OUT]
                         LAND_RESTART_IN_BASENAME, 'ZtS', 'XY', step=1    ) ! [IN]
       call FILEIO_read( Ze0(:,:),                                        & ! [OUT]
                         LAND_RESTART_IN_BASENAME, 'Ze0', 'XY', step=1    ) ! [IN]
       call FILEIO_read( ZeR(:,:),                                        & ! [OUT]
                         LAND_RESTART_IN_BASENAME, 'ZeR', 'XY', step=1    ) ! [IN]
       call FILEIO_read( ZeS(:,:),                                        & ! [OUT]
                         LAND_RESTART_IN_BASENAME, 'ZeS', 'XY', step=1    ) ! [IN]

!       call FILEIO_read( SoilT(:,:,:),                                    & ! [OUT]
!                         LAND_RESTART_IN_BASENAME, 'SoilT', 'ZXY', step=1 ) ! [IN]
!       call FILEIO_read( SoilW(:,:,:),                                    & ! [OUT]
!                         LAND_RESTART_IN_BASENAME, 'SoilW', 'ZXY', step=1 ) ! [IN]
!       call FILEIO_read( SoilI(:,:,:),                                    & ! [OUT]
!                         LAND_RESTART_IN_BASENAME, 'SoilI', 'ZXY', step=1 ) ! [IN]

       call LAND_vars_fillhalo

       call LAND_vars_total

    else

       if( IO_L ) write(IO_FID_LOG,*) '*** restart file for land is not specified.'

!       TG   (:,:) = CONST_UNDEF
!       QvEfc(:,:) = CONST_UNDEF
!       EMIT (:,:) = CONST_UNDEF
!       ALB  (:,:) = CONST_UNDEF
!       TCS  (:,:) = CONST_UNDEF
!       HCS  (:,:) = CONST_UNDEF
!       DZg  (:,:) = CONST_UNDEF
!       Z00  (:,:) = CONST_UNDEF
!       Z0R  (:,:) = CONST_UNDEF
!       Z0S  (:,:) = CONST_UNDEF
!       Zt0  (:,:) = CONST_UNDEF
!       ZtR  (:,:) = CONST_UNDEF
!       ZtS  (:,:) = CONST_UNDEF
!       Ze0  (:,:) = CONST_UNDEF
!       ZeR  (:,:) = CONST_UNDEF
!       ZeS  (:,:) = CONST_UNDEF

       TG   (:,:) = 300.0_RP
       QvEfc(:,:) = 1.0_RP
       EMIT (:,:) = 0.98_RP
       ALB  (:,:) = 0.33_RP
       TCS  (:,:) = 1.0_RP
       HCS  (:,:) = 2.2E+6_RP
       DZg  (:,:) = 1.0_RP
       Z00  (:,:) = 0.0_RP
       Z0R  (:,:) = 0.018_RP
       Z0S  (:,:) = 0.11_RP
       Zt0  (:,:) = 1.4E-5_RP
       ZtR  (:,:) = 0.0_RP
       ZtS  (:,:) = 0.4_RP
       Ze0  (:,:) = 1.3E-4_RP
       ZeR  (:,:) = 0.0_RP
       ZeS  (:,:) = 0.62_RP

!       SoilT(:,:,:) = CONST_UNDEF
!       SoilW(:,:,:) = CONST_UNDEF
!       SoilI(:,:,:) = CONST_UNDEF

    endif

    call TIME_rapend  ('FILE I NetCDF')

    return
  end subroutine LAND_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Write land restart
  subroutine LAND_vars_restart_write
    use mod_time, only: &
       NOWSEC => TIME_NOWDAYSEC
    use mod_fileio, only: &
       FILEIO_write
    implicit none

    character(len=IO_FILECHR) :: bname

    integer :: n
    !---------------------------------------------------------------------------

    call TIME_rapstart('FILE O NetCDF')

    if ( LAND_RESTART_OUT_BASENAME /= '' ) then

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Output restart file (land) ***'

       bname = ''
       write(bname(1:15), '(F15.3)') NOWSEC
       do n = 1, 15
          if ( bname(n:n) == ' ' ) bname(n:n) = '0'
       enddo
       write(bname,'(A,A,A)') trim(LAND_RESTART_OUT_BASENAME), '_', trim(bname)

       call FILEIO_write( TG(:,:), bname, LAND_RESTART_OUT_TITLE,                                              & ! [IN]
                          LP_NAME(I_TG),    LP_DESC(I_TG),    LP_UNIT(I_TG),    'XY', LAND_RESTART_OUT_DTYPE   ) ! [IN]
       call FILEIO_write( QvEfc(:,:), bname, LAND_RESTART_OUT_TITLE,                                           & ! [IN]
                          LP_NAME(I_QvEfc), LP_DESC(I_QvEfc), LP_UNIT(I_QvEfc), 'XY', LAND_RESTART_OUT_DTYPE   ) ! [IN]
       call FILEIO_write( EMIT(:,:), bname, LAND_RESTART_OUT_TITLE,                                            & ! [IN]
                          LP_NAME(I_EMIT),  LP_DESC(I_EMIT),  LP_UNIT(I_EMIT),  'XY', LAND_RESTART_OUT_DTYPE   ) ! [IN]
       call FILEIO_write( ALB(:,:), bname, LAND_RESTART_OUT_TITLE,                                             & ! [IN]
                          LP_NAME(I_ALB),   LP_DESC(I_ALB),   LP_UNIT(I_ALB),   'XY', LAND_RESTART_OUT_DTYPE   ) ! [IN]
       call FILEIO_write( TCS(:,:), bname, LAND_RESTART_OUT_TITLE,                                             & ! [IN]
                          LP_NAME(I_TCS),   LP_DESC(I_TCS),   LP_UNIT(I_TCS),   'XY', LAND_RESTART_OUT_DTYPE   ) ! [IN]
       call FILEIO_write( HCS(:,:), bname, LAND_RESTART_OUT_TITLE,                                             & ! [IN]
                          LP_NAME(I_HCS),   LP_DESC(I_HCS),   LP_UNIT(I_HCS),   'XY', LAND_RESTART_OUT_DTYPE   ) ! [IN]
       call FILEIO_write( DZg(:,:), bname, LAND_RESTART_OUT_TITLE,                                             & ! [IN]
                          LP_NAME(I_DZg),   LP_DESC(I_DZg),   LP_UNIT(I_DZg),   'XY', LAND_RESTART_OUT_DTYPE   ) ! [IN]
       call FILEIO_write( Z00(:,:), bname, LAND_RESTART_OUT_TITLE,                                             & ! [IN]
                          LP_NAME(I_Z00),   LP_DESC(I_Z00),   LP_UNIT(I_Z00),   'XY', LAND_RESTART_OUT_DTYPE   ) ! [IN]
       call FILEIO_write( Z0R(:,:), bname, LAND_RESTART_OUT_TITLE,                                             & ! [IN]
                          LP_NAME(I_Z0R),   LP_DESC(I_Z0R),   LP_UNIT(I_Z0R),   'XY', LAND_RESTART_OUT_DTYPE   ) ! [IN]
       call FILEIO_write( Z0S(:,:), bname, LAND_RESTART_OUT_TITLE,                                             & ! [IN]
                          LP_NAME(I_Z0S),   LP_DESC(I_Z0S),   LP_UNIT(I_Z0S),   'XY', LAND_RESTART_OUT_DTYPE   ) ! [IN]
       call FILEIO_write( Zt0(:,:), bname, LAND_RESTART_OUT_TITLE,                                             & ! [IN]
                          LP_NAME(I_Zt0),   LP_DESC(I_Zt0),   LP_UNIT(I_Zt0),   'XY', LAND_RESTART_OUT_DTYPE   ) ! [IN]
       call FILEIO_write( ZtR(:,:), bname, LAND_RESTART_OUT_TITLE,                                             & ! [IN]
                          LP_NAME(I_ZtR),   LP_DESC(I_ZtR),   LP_UNIT(I_ZtR),   'XY', LAND_RESTART_OUT_DTYPE   ) ! [IN]
       call FILEIO_write( ZtS(:,:), bname, LAND_RESTART_OUT_TITLE,                                             & ! [IN]
                          LP_NAME(I_ZtS),   LP_DESC(I_ZtS),   LP_UNIT(I_ZtS),   'XY', LAND_RESTART_OUT_DTYPE   ) ! [IN]
       call FILEIO_write( Ze0(:,:), bname, LAND_RESTART_OUT_TITLE,                                             & ! [IN]
                          LP_NAME(I_Ze0),   LP_DESC(I_Ze0),   LP_UNIT(I_Ze0),   'XY', LAND_RESTART_OUT_DTYPE   ) ! [IN]
       call FILEIO_write( ZeR(:,:), bname, LAND_RESTART_OUT_TITLE,                                             & ! [IN]
                          LP_NAME(I_ZeR),   LP_DESC(I_ZeR),   LP_UNIT(I_ZeR),   'XY', LAND_RESTART_OUT_DTYPE   ) ! [IN]
       call FILEIO_write( ZeS(:,:), bname, LAND_RESTART_OUT_TITLE,                                             & ! [IN]
                          LP_NAME(I_ZeS),   LP_DESC(I_ZeS),   LP_UNIT(I_ZeS),   'XY', LAND_RESTART_OUT_DTYPE   ) ! [IN]

!       call FILEIO_write( SoilT(:,:,:), bname, LAND_RESTART_OUT_TITLE,                                         & ! [IN]
!                          LP_NAME(I_SoilT), LP_DESC(I_SoilT), LP_UNIT(I_SoilT), 'ZXY', LAND_RESTART_OUT_DTYPE  ) ! [IN]
!       call FILEIO_write( SoilW(:,:,:), bname, LAND_RESTART_OUT_TITLE,                                         & ! [IN]
!                          LP_NAME(I_SoilW), LP_DESC(I_SoilW), LP_UNIT(I_SoilW), 'ZXY', LAND_RESTART_OUT_DTYPE  ) ! [IN]
!       call FILEIO_write( SoilI(:,:,:), bname, LAND_RESTART_OUT_TITLE,                                         & ! [IN]
!                          LP_NAME(I_SoilI), LP_DESC(I_SoilI), LP_UNIT(I_SoilI), 'ZXY', LAND_RESTART_OUT_DTYPE  ) ! [IN]

    endif

    call TIME_rapend  ('FILE O NetCDF')

    call LAND_vars_total

    return
  end subroutine LAND_vars_restart_write

  !-----------------------------------------------------------------------------
  !> History output set for land variables
  subroutine LAND_vars_history
    use mod_misc, only: &
       MISC_valcheck
    use mod_time, only: &
       TIME_DTSEC_LAND
    use mod_history, only: &
       HIST_in
    implicit none
    !---------------------------------------------------------------------------

    if ( LAND_VARS_CHECKRANGE ) then
       call MISC_valcheck( TG   (:,:),      0.0_RP, 1000.0_RP,  LP_NAME(I_TG)    )
       call MISC_valcheck( QvEfc(:,:),      0.0_RP,    2.0_RP,  LP_NAME(I_QvEfc) )
       call MISC_valcheck( EMIT (:,:),      0.0_RP,    2.0_RP,  LP_NAME(I_EMIT)  )
       call MISC_valcheck( ALB  (:,:),      0.0_RP,    2.0_RP,  LP_NAME(I_ALB)   )
       call MISC_valcheck( TCS  (:,:),      0.0_RP, 1000.0_RP,  LP_NAME(I_TCS)   )
       call MISC_valcheck( HCS  (:,:),      0.0_RP, 1.0E+10_RP, LP_NAME(I_HCS)   )
       call MISC_valcheck( DZg  (:,:),      0.0_RP, 1000.0_RP,  LP_NAME(I_DZg)   )
       call MISC_valcheck( Z00  (:,:),      0.0_RP, 1000.0_RP,  LP_NAME(I_Z00)   )
       call MISC_valcheck( Z0R  (:,:),      0.0_RP, 1000.0_RP,  LP_NAME(I_Z0R)   )
       call MISC_valcheck( Z0S  (:,:),      0.0_RP, 1000.0_RP,  LP_NAME(I_Z0S)   )
       call MISC_valcheck( Zt0  (:,:),      0.0_RP, 1000.0_RP,  LP_NAME(I_Zt0)   )
       call MISC_valcheck( ZtR  (:,:),      0.0_RP, 1000.0_RP,  LP_NAME(I_ZtR)   )
       call MISC_valcheck( ZtS  (:,:),      0.0_RP, 1000.0_RP,  LP_NAME(I_ZtS)   )
       call MISC_valcheck( Ze0  (:,:),      0.0_RP, 1000.0_RP,  LP_NAME(I_Ze0)   )
       call MISC_valcheck( ZeR  (:,:),      0.0_RP, 1000.0_RP,  LP_NAME(I_ZeR)   )
       call MISC_valcheck( ZeS  (:,:),      0.0_RP, 1000.0_RP,  LP_NAME(I_ZeS)   )
!       call MISC_valcheck( SoilT(:,:,:),    0.0_RP, 1000.0_RP,  LP_NAME(I_SoilT) )
!       call MISC_valcheck( SoilW(:,:,:),    0.0_RP,    2.0_RP,  LP_NAME(I_SoilW) )
!       call MISC_valcheck( SoilI(:,:,:),    0.0_RP,    2.0_RP,  LP_NAME(I_SoilI) )
    endif

    call HIST_in( TG   (:,:),   'L_TG',    LP_DESC(I_TG),    LP_UNIT(I_TG),    TIME_DTSEC_LAND )
    call HIST_in( QvEfc(:,:),   'L_QvEfc', LP_DESC(I_QvEfc), LP_UNIT(I_QvEfc), TIME_DTSEC_LAND )
    call HIST_in( EMIT (:,:),   'L_EMIT',  LP_DESC(I_EMIT),  LP_UNIT(I_EMIT),  TIME_DTSEC_LAND )
    call HIST_in( ALB  (:,:),   'L_ALB',   LP_DESC(I_ALB),   LP_UNIT(I_ALB),   TIME_DTSEC_LAND )
    call HIST_in( TCS  (:,:),   'L_TCS',   LP_DESC(I_TCS),   LP_UNIT(I_TCS),   TIME_DTSEC_LAND )
    call HIST_in( HCS  (:,:),   'L_HCS',   LP_DESC(I_HCS),   LP_UNIT(I_HCS),   TIME_DTSEC_LAND )
    call HIST_in( DZg  (:,:),   'L_DZg',   LP_DESC(I_DZg),   LP_UNIT(I_DZg),   TIME_DTSEC_LAND )
    call HIST_in( Z00  (:,:),   'L_Z00',   LP_DESC(I_Z00),   LP_UNIT(I_Z00),   TIME_DTSEC_LAND )
    call HIST_in( Z0R  (:,:),   'L_Z0R',   LP_DESC(I_Z0R),   LP_UNIT(I_Z0R),   TIME_DTSEC_LAND )
    call HIST_in( Z0S  (:,:),   'L_Z0S',   LP_DESC(I_Z0S),   LP_UNIT(I_Z0S),   TIME_DTSEC_LAND )
    call HIST_in( Zt0  (:,:),   'L_Zt0',   LP_DESC(I_Zt0),   LP_UNIT(I_Zt0),   TIME_DTSEC_LAND )
    call HIST_in( ZtR  (:,:),   'L_ZtR',   LP_DESC(I_ZtR),   LP_UNIT(I_ZtR),   TIME_DTSEC_LAND )
    call HIST_in( ZtS  (:,:),   'L_ZtS',   LP_DESC(I_ZtS),   LP_UNIT(I_ZtS),   TIME_DTSEC_LAND )
    call HIST_in( Ze0  (:,:),   'L_Ze0',   LP_DESC(I_Ze0),   LP_UNIT(I_Ze0),   TIME_DTSEC_LAND )
    call HIST_in( ZeR  (:,:),   'L_ZeR',   LP_DESC(I_ZeR),   LP_UNIT(I_ZeR),   TIME_DTSEC_LAND )
    call HIST_in( ZeS  (:,:),   'L_ZeS',   LP_DESC(I_ZeS),   LP_UNIT(I_ZeS),   TIME_DTSEC_LAND )
!    call HIST_in( SoilT(:,:,:), 'L_SoilT', LP_DESC(I_SoilT), LP_UNIT(I_SoilT), TIME_DTSEC_LAND )
!    call HIST_in( SoilW(:,:,:), 'L_SoilW', LP_DESC(I_SoilW), LP_UNIT(I_SoilW), TIME_DTSEC_LAND )
!    call HIST_in( SoilI(:,:,:), 'L_SoilI', LP_DESC(I_SoilI), LP_UNIT(I_SoilI), TIME_DTSEC_LAND )

    return
  end subroutine LAND_vars_history

  !-----------------------------------------------------------------------------
  !> Budget monitor for land
  subroutine LAND_vars_total
    use mod_comm, only: &
       COMM_total_doreport, &
       COMM_total
    implicit none

    real(RP) :: total ! dummy
    !---------------------------------------------------------------------------

!    if ( COMM_total_doreport ) then
!
!       call COMM_total( total, SoilT(:,:,:), LP_NAME(I_SoilT) )
!       call COMM_total( total, SoilW(:,:,:), LP_NAME(I_SoilW) )
!       call COMM_total( total, SoilI(:,:,:), LP_NAME(I_SoilI) )
!
!    endif

    return
  end subroutine LAND_vars_total

end module mod_land_vars
