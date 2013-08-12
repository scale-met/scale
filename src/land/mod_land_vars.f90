!-------------------------------------------------------------------------------
!> module LAND Variables
!!
!! @par Description
!!          Container for land variables
!!
!! @author Team SCALE
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
  include 'inc_land.h'

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
  real(RP), public, save :: SFLX_GH  (IA,JA) ! ground heat flux (upward positive) [W/m2]
  real(RP), public, save :: SFLX_PREC(IA,JA) ! precipitation flux [kg/m2/s]
  real(RP), public, save :: SFLX_QV  (IA,JA) ! moisture flux [kg/m2/s]

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

  real(RP), public, save :: SkinT(IA,JA) ! Ground Skin Temperature [K]
  real(RP), public, save :: SkinW(IA,JA) ! Ground Skin Water       [kg/m2]
  real(RP), public, save :: SnowQ(IA,JA) ! Ground Snow amount      [kg/m2]
  real(RP), public, save :: SnowT(IA,JA) ! Ground Snow Temperature [K]

  real(RP), public, save :: SoilT(KA_soil,IA,JA) ! Soil temperature             [K]
  real(RP), public, save :: SoilW(KA_soil,IA,JA) ! Soil moisture (liquid water) [m3/m3]
  real(RP), public, save :: SoilI(KA_soil,IA,JA) ! Soil ice      (ice    water) [m3/m3]

  integer,                    public, save :: I_SkinT = 1
  integer,                    public, save :: I_SkinW = 1
  integer,                    public, save :: I_SnowQ = 1
  integer,                    public, save :: I_SnowT = 1
  integer,                    public, save :: I_SoilT = 1
  integer,                    public, save :: I_SoilW = 1
  integer,                    public, save :: I_SoilI = 1
  character(len=File_HSHORT), public, save :: LP_NAME(7) !< name  of the land variables
  character(len=File_HMID),   public, save :: LP_DESC(7) !< desc. of the land variables
  character(len=File_HSHORT), public, save :: LP_UNIT(7) !< unit  of the land variables

  data LP_NAME / 'SkinT', &
                 'SkinW', &
                 'SnowQ', &
                 'SnowT', &
                 'SoilT', &
                 'SoilW', &
                 'SoilI'  /
  data LP_DESC / 'ground skin temp.',  &
                 'ground skin water',  &
                 'ground snow amount', &
                 'ground snow temp.',  &
                 'soil temperature',   &
                 'soil moisture',      &
                 'soil ice'            /
  data LP_UNIT / 'K',     &
                 'kg/m2', &
                 'kg/m2', &
                 'K',     &
                 'K',     &
                 'm3/m3', &
                 'm3/m3'  /

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
    do ip = 1, 7
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

    TG   (:,:) = 285.0_RP
    QvEfc(:,:) = 1.0_RP
    EMIT (:,:) = 0.98_RP
    ALB  (:,:) = 0.33_RP
    TCS  (:,:) = 1.0_RP
    HCS  (:,:) = 2.2E6_RP
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
    call COMM_vars8( SkinT(:,:), 1 )
    call COMM_vars8( SkinW(:,:), 2 )

    call COMM_vars8( SnowQ(:,:), 3 )
    call COMM_vars8( SnowT(:,:), 4 )

    call COMM_vars8( SoilT(:,:,:), 5 )
    call COMM_vars8( SoilW(:,:,:), 6 )
    call COMM_vars8( SoilI(:,:,:), 7 )

    call COMM_wait ( SkinT(:,:), 1 )
    call COMM_wait ( SkinW(:,:), 2 )

    call COMM_wait ( SnowQ(:,:), 3 )
    call COMM_wait ( SnowT(:,:), 4 )

    call COMM_wait ( SoilT(:,:,:), 5 )
    call COMM_wait ( SoilW(:,:,:), 6 )
    call COMM_wait ( SoilI(:,:,:), 7 )

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

       call FILEIO_read( SkinT(:,:),                                     & ! [OUT]
                         LAND_RESTART_IN_BASENAME, 'SkinT', 'XY', step=1 ) ! [IN]
       call FILEIO_read( SkinW(:,:),                                     & ! [OUT]
                         LAND_RESTART_IN_BASENAME, 'SkinW', 'XY', step=1 ) ! [IN]
       call FILEIO_read( SnowQ(:,:),                                     & ! [OUT]
                         LAND_RESTART_IN_BASENAME, 'SnowQ', 'XY', step=1 ) ! [IN]
       call FILEIO_read( SnowT(:,:),                                     & ! [OUT]
                         LAND_RESTART_IN_BASENAME, 'SnowT', 'XY', step=1 ) ! [IN]

       call FILEIO_read( SoilT(:,:,:),                                    & ! [OUT]
                         LAND_RESTART_IN_BASENAME, 'SoilT', 'ZXY', step=1 ) ! [IN]
       call FILEIO_read( SoilW(:,:,:),                                    & ! [OUT]
                         LAND_RESTART_IN_BASENAME, 'SoilW', 'ZXY', step=1 ) ! [IN]
       call FILEIO_read( SoilI(:,:,:),                                    & ! [OUT]
                         LAND_RESTART_IN_BASENAME, 'SoilI', 'ZXY', step=1 ) ! [IN]

       call LAND_vars_fillhalo

       call LAND_vars_total
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** restart file for land is not specified.'

       SkinT(:,:) = CONST_UNDEF
       SkinW(:,:) = CONST_UNDEF
       SnowQ(:,:) = CONST_UNDEF
       SnowT(:,:) = CONST_UNDEF

       SoilT(:,:,:) = CONST_UNDEF
       SoilW(:,:,:) = CONST_UNDEF
       SoilI(:,:,:) = CONST_UNDEF
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

       call FILEIO_write( SkinT(:,:), bname,                        LAND_RESTART_OUT_TITLE, & ! [IN]
                          LP_NAME(1), LP_DESC(1), LP_UNIT(1), 'XY', LAND_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( SkinW(:,:), bname,                        LAND_RESTART_OUT_TITLE, & ! [IN]
                          LP_NAME(2), LP_DESC(2), LP_UNIT(2), 'XY', LAND_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( SnowQ(:,:), bname,                        LAND_RESTART_OUT_TITLE, & ! [IN]
                          LP_NAME(3), LP_DESC(3), LP_UNIT(3), 'XY', LAND_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( SnowT(:,:), bname,                        LAND_RESTART_OUT_TITLE, & ! [IN]
                          LP_NAME(4), LP_DESC(4), LP_UNIT(4), 'XY', LAND_RESTART_OUT_DTYPE  ) ! [IN]

       call FILEIO_write( SoilT(:,:,:), bname,                       LAND_RESTART_OUT_TITLE, & ! [IN]
                          LP_NAME(5), LP_DESC(5), LP_UNIT(5), 'ZXY', LAND_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( SoilW(:,:,:), bname,                       LAND_RESTART_OUT_TITLE, & ! [IN]
                          LP_NAME(6), LP_DESC(6), LP_UNIT(6), 'ZXY', LAND_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( SoilI(:,:,:), bname,                       LAND_RESTART_OUT_TITLE, & ! [IN]
                          LP_NAME(7), LP_DESC(7), LP_UNIT(7), 'ZXY', LAND_RESTART_OUT_DTYPE  ) ! [IN]

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
       call MISC_valcheck( SkinT(:,:),      0.0_RP, 1000.0_RP, LP_NAME(I_SkinT) )
       call MISC_valcheck( SkinW(:,:),      0.0_RP, 1000.0_RP, LP_NAME(I_SkinW) )
       call MISC_valcheck( SnowQ(:,:),      0.0_RP, 1000.0_RP, LP_NAME(I_SnowQ) )
       call MISC_valcheck( SnowT(:,:),      0.0_RP, 1000.0_RP, LP_NAME(I_SnowT) )
       call MISC_valcheck( SoilT(:,:,:),    0.0_RP, 1000.0_RP, LP_NAME(I_SoilT) )
       call MISC_valcheck( SoilW(:,:,:),    0.0_RP,    2.0_RP, LP_NAME(I_SoilW) )
       call MISC_valcheck( SoilI(:,:,:),    0.0_RP,    2.0_RP, LP_NAME(I_SoilI) )
    endif

    call HIST_in( SkinT(:,:),   'L_SkinT', LP_DESC(I_SkinT), LP_UNIT(I_SkinT), TIME_DTSEC_LAND )
    call HIST_in( SkinW(:,:),   'L_SkinW', LP_DESC(I_SkinW), LP_UNIT(I_SkinW), TIME_DTSEC_LAND )
    call HIST_in( SnowQ(:,:),   'L_SnowQ', LP_DESC(I_SnowQ), LP_UNIT(I_SnowQ), TIME_DTSEC_LAND )
    call HIST_in( SnowT(:,:),   'L_SnowT', LP_DESC(I_SnowT), LP_UNIT(I_SnowT), TIME_DTSEC_LAND )
    call HIST_in( SoilT(:,:,:), 'L_SoilT', LP_DESC(I_SoilT), LP_UNIT(I_SoilT), TIME_DTSEC_LAND )
    call HIST_in( SoilW(:,:,:), 'L_SoilW', LP_DESC(I_SoilW), LP_UNIT(I_SoilW), TIME_DTSEC_LAND )
    call HIST_in( SoilI(:,:,:), 'L_SoilI', LP_DESC(I_SoilI), LP_UNIT(I_SoilI), TIME_DTSEC_LAND )

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
!       call COMM_total( total, SkinT(:,:),   LP_NAME(I_SkinT) )
!       call COMM_total( total, SkinW(:,:),   LP_NAME(I_SkinW) )
!       call COMM_total( total, SnowQ(:,:),   LP_NAME(I_SnowQ) )
!       call COMM_total( total, SnowT(:,:),   LP_NAME(I_SnowT) )
!       call COMM_total( total, SoilT(:,:,:), LP_NAME(I_SoilT) )
!       call COMM_total( total, SoilW(:,:,:), LP_NAME(I_SoilW) )
!       call COMM_total( total, SoilI(:,:,:), LP_NAME(I_SoilI) )
!
!    endif

    return
  end subroutine LAND_vars_total

end module mod_land_vars
