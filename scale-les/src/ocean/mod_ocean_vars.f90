!-------------------------------------------------------------------------------
!> module OCEAN Variables
!!
!! @par Description
!!          Container for ocean variables
!!
!! @author Team SCALE
!! @li      2013-08-31 (T.Yamaura)  [new]
!<
!-------------------------------------------------------------------------------
module mod_ocean_vars
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_debug
  use scale_grid_index

  use scale_const, only: &
     I_SW  => CONST_I_SW, &
     I_LW  => CONST_I_LW
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: OCEAN_vars_setup
  public :: OCEAN_vars_restart_read
  public :: OCEAN_vars_restart_write
  public :: OCEAN_vars_history
  public :: OCEAN_vars_total
  public :: OCEAN_vars_external_in

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical,               public :: OCEAN_RESTART_OUTPUT       = .false.         !< output restart file?

  character(len=H_LONG), public :: OCEAN_RESTART_IN_BASENAME  = ''              !< basename of the input file
  character(len=H_LONG), public :: OCEAN_RESTART_OUT_BASENAME = ''              !< basename of the output file
  character(len=H_MID),  public :: OCEAN_RESTART_OUT_TITLE    = 'OCEAN restart' !< title    of the output file
  character(len=H_MID),  public :: OCEAN_RESTART_OUT_DTYPE    = 'DEFAULT'       !< REAL4 or REAL8

  ! prognostic variables
  real(RP), public, allocatable :: OCEAN_TEMP      (:,:)   !< temperature at uppermost ocean layer (SST)  [K]

  ! tendency variables
  real(RP), public, allocatable :: OCEAN_TEMP_t    (:,:)   !< tendency of OCEAN_TEMP

  ! for restart
  real(RP), public, allocatable :: OCEAN_SFC_TEMP  (:,:)   !< ocean surface skin temperature              [K]
  real(RP), public, allocatable :: OCEAN_SFC_albedo(:,:,:) !< ocean surface albedo                        [0-1]
  real(RP), public, allocatable :: OCEAN_SFC_Z0M   (:,:)   !< ocean surface roughness length for momentum [m]
  real(RP), public, allocatable :: OCEAN_SFC_Z0H   (:,:)   !< ocean surface roughness length for heat [m]
  real(RP), public, allocatable :: OCEAN_SFC_Z0E   (:,:)   !< ocean surface roughness length for vapor [m]

  real(RP), public, allocatable :: OCEAN_SFLX_heat (:,:)   !< ocean surface heat flux [J/m2/s]
  real(RP), public, allocatable :: OCEAN_SFLX_prec (:,:)   !< ocean surface precipitation flux [kg/m2/s]
  real(RP), public, allocatable :: OCEAN_SFLX_evap (:,:)   !< ocean surface water vapor flux [kg/m2/s]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical,                private :: OCEAN_VARS_CHECKRANGE      = .false.

  integer,                private, parameter :: VMAX = 10      !< number of the variables
  integer,                private, parameter :: I_TEMP          =  1
  integer,                private, parameter :: I_SFC_TEMP      =  2
  integer,                private, parameter :: I_SFC_albedo_LW =  3
  integer,                private, parameter :: I_SFC_albedo_SW =  4
  integer,                private, parameter :: I_SFC_Z0M       =  5
  integer,                private, parameter :: I_SFC_Z0H       =  6
  integer,                private, parameter :: I_SFC_Z0E       =  7
  integer,                private, parameter :: I_SFLX_heat     =  8
  integer,                private, parameter :: I_SFLX_prec     =  9
  integer,                private, parameter :: I_SFLX_evap     = 10

  character(len=H_SHORT), private            :: VAR_NAME(VMAX) !< name  of the variables
  character(len=H_MID),   private            :: VAR_DESC(VMAX) !< desc. of the variables
  character(len=H_SHORT), private            :: VAR_UNIT(VMAX) !< unit  of the variables

  data VAR_NAME / 'OCEAN_TEMP',      &
                  'OCEAN_SFC_TEMP',  &
                  'OCEAN_ALB_LW',    &
                  'OCEAN_ALB_SW',    &
                  'OCEAN_SFC_Z0M',   &
                  'OCEAN_SFC_Z0H',   &
                  'OCEAN_SFC_Z0E',   &
                  'OCEAN_SFLX_heat', &
                  'OCEAN_SFLX_prec', &
                  'OCEAN_SFLX_evap'  /
  data VAR_DESC / 'temperature at uppermost ocean layer (SST)', &
                  'ocean surface skin temperature',             &
                  'ocean surface albedo (longwave)',            &
                  'ocean surface albedo (shortwave)',           &
                  'ocean surface roughness length (momentum)',  &
                  'ocean surface roughness length (heat)',      &
                  'ocean surface roughness length (vapor)',     &
                  'ocean surface heat flux',                    &
                  'ocean surface precipitation flux',           &
                  'ocean surface water vapor flux'              /
  data VAR_UNIT / 'K',       &
                  'K',       &
                  '0-1',     &
                  '0-1',     &
                  'm',       &
                  'm',       &
                  'm',       &
                  'J/m2/s',  &
                  'kg/m2/s', &
                  'kg/m2/s'  /

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine OCEAN_vars_setup
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    implicit none

    NAMELIST / PARAM_OCEAN_VARS /  &
       OCEAN_RESTART_IN_BASENAME,  &
       OCEAN_RESTART_OUTPUT,       &
       OCEAN_RESTART_OUT_BASENAME, &
       OCEAN_RESTART_OUT_TITLE,    &
       OCEAN_VARS_CHECKRANGE

    integer :: ierr
    integer :: iv
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[VARS] / Categ[OCEAN] / Origin[SCALE-LES]'

    allocate( OCEAN_TEMP  (IA,JA) )
    allocate( OCEAN_TEMP_t(IA,JA) )
    OCEAN_TEMP  (:,:) = UNDEF
    OCEAN_TEMP_t(:,:) = UNDEF

    allocate( OCEAN_SFC_TEMP  (IA,JA)   )
    allocate( OCEAN_SFC_albedo(IA,JA,2) )
    allocate( OCEAN_SFC_Z0M   (IA,JA)   )
    allocate( OCEAN_SFC_Z0H   (IA,JA)   )
    allocate( OCEAN_SFC_Z0E   (IA,JA)   )
    OCEAN_SFC_TEMP  (:,:)   = UNDEF
    OCEAN_SFC_albedo(:,:,:) = UNDEF
    OCEAN_SFC_Z0M   (:,:)   = UNDEF
    OCEAN_SFC_Z0H   (:,:)   = UNDEF
    OCEAN_SFC_Z0E   (:,:)   = UNDEF

    allocate( OCEAN_SFLX_heat(IA,JA) )
    allocate( OCEAN_SFLX_prec(IA,JA) )
    allocate( OCEAN_SFLX_evap(IA,JA) )
    OCEAN_SFLX_heat(:,:) = UNDEF
    OCEAN_SFLX_prec(:,:) = UNDEF
    OCEAN_SFLX_evap(:,:) = UNDEF

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCEAN_VARS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_OCEAN_VARS. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_OCEAN_VARS)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** List of prognostic variables (OCEAN) ***'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A15,A,A32,3(A))') &
               '***       |','VARNAME        ','|', 'DESCRIPTION                     ','[', 'UNIT            ',']'
    do iv = 1, VMAX
       if( IO_L ) write(IO_FID_LOG,'(1x,A,i3,A,A15,A,A32,3(A))') &
                  '*** NO.',iv,'|',VAR_NAME(iv),'|',VAR_DESC(iv),'[',VAR_UNIT(iv),']'

    enddo

    if( IO_L ) write(IO_FID_LOG,*)
    if ( OCEAN_RESTART_IN_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart input?  : ', trim(OCEAN_RESTART_IN_BASENAME)
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart input?  : NO'
    endif
    if (       OCEAN_RESTART_OUTPUT             &
         .AND. OCEAN_RESTART_OUT_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart output? : ', trim(OCEAN_RESTART_OUT_BASENAME)
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart output? : NO'
       OCEAN_RESTART_OUTPUT = .false.
    endif

    return
  end subroutine OCEAN_vars_setup

  !-----------------------------------------------------------------------------
  !> Read ocean restart
  subroutine OCEAN_vars_restart_read
    use scale_fileio, only: &
       FILEIO_read
    use mod_ocean_admin, only: &
       OCEAN_sw
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input restart file (OCEAN) ***'

    if ( OCEAN_sw .and. OCEAN_RESTART_IN_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(OCEAN_RESTART_IN_BASENAME)

       call FILEIO_read( OCEAN_TEMP(:,:),                                                   & ! [OUT]
                         OCEAN_RESTART_IN_BASENAME, VAR_NAME(I_TEMP),          'XY', step=1 ) ! [IN]
       call FILEIO_read( OCEAN_SFC_TEMP(:,:),                                               & ! [OUT]
                         OCEAN_RESTART_IN_BASENAME, VAR_NAME(I_SFC_TEMP),      'XY', step=1 ) ! [IN]
       call FILEIO_read( OCEAN_SFC_albedo(:,:,I_LW),                                        & ! [OUT]
                         OCEAN_RESTART_IN_BASENAME, VAR_NAME(I_SFC_albedo_LW), 'XY', step=1 ) ! [IN]
       call FILEIO_read( OCEAN_SFC_albedo(:,:,I_SW),                                        & ! [OUT]
                         OCEAN_RESTART_IN_BASENAME, VAR_NAME(I_SFC_albedo_SW), 'XY', step=1 ) ! [IN]
       call FILEIO_read( OCEAN_SFC_Z0M(:,:),                                                & ! [OUT]
                         OCEAN_RESTART_IN_BASENAME, VAR_NAME(I_SFC_Z0M),       'XY', step=1 ) ! [IN]
       call FILEIO_read( OCEAN_SFC_Z0H(:,:),                                                & ! [OUT]
                         OCEAN_RESTART_IN_BASENAME, VAR_NAME(I_SFC_Z0H),       'XY', step=1 ) ! [IN]
       call FILEIO_read( OCEAN_SFC_Z0E(:,:),                                                & ! [OUT]
                         OCEAN_RESTART_IN_BASENAME, VAR_NAME(I_SFC_Z0E),       'XY', step=1 ) ! [IN]
       call FILEIO_read( OCEAN_SFLX_heat(:,:),                                              & ! [OUT]
                         OCEAN_RESTART_IN_BASENAME, VAR_NAME(I_SFLX_heat),     'XY', step=1 ) ! [IN]
       call FILEIO_read( OCEAN_SFLX_prec(:,:),                                              & ! [OUT]
                         OCEAN_RESTART_IN_BASENAME, VAR_NAME(I_SFLX_prec),     'XY', step=1 ) ! [IN]
       call FILEIO_read( OCEAN_SFLX_evap(:,:),                                              & ! [OUT]
                         OCEAN_RESTART_IN_BASENAME, VAR_NAME(I_SFLX_evap),     'XY', step=1 ) ! [IN]

       call OCEAN_vars_total
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** restart file for ocean is not specified.'
    endif

    return
  end subroutine OCEAN_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Write ocean restart
  subroutine OCEAN_vars_restart_write
    use scale_time, only: &
       TIME_gettimelabel
    use scale_fileio, only: &
       FILEIO_write
    use mod_ocean_admin, only: &
       OCEAN_sw
    use mod_cpl_admin, only: &
       CPL_sw_AtmOcn
    use mod_cpl_vars, only: &
       CPL_getOcn
    implicit none

    character(len=15)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    if ( OCEAN_sw .and. OCEAN_RESTART_OUT_BASENAME /= '' ) then

       if ( CPL_sw_AtmOcn ) then
          call CPL_getOcn( OCEAN_SFC_TEMP  (:,:),   & ! [OUT]
                           OCEAN_SFC_albedo(:,:,:), & ! [OUT]
                           OCEAN_SFC_Z0M   (:,:),   & ! [OUT]
                           OCEAN_SFC_Z0H   (:,:),   & ! [OUT]
                           OCEAN_SFC_Z0E   (:,:),   & ! [OUT]
                           OCEAN_SFLX_heat (:,:),   & ! [OUT]
                           OCEAN_SFLX_prec (:,:),   & ! [OUT]
                           OCEAN_SFLX_evap (:,:)    ) ! [OUT]
       endif

       call TIME_gettimelabel( timelabel )
       write(basename,'(A,A,A)') trim(OCEAN_RESTART_OUT_BASENAME), '_', trim(timelabel)

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Output restart file (OCEAN) ***'
       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(basename)

       call OCEAN_vars_total

       call FILEIO_write( OCEAN_TEMP(:,:),            basename,                  OCEAN_RESTART_OUT_TITLE,   & ! [IN]
                          VAR_NAME(I_TEMP),           VAR_DESC(I_TEMP),          VAR_UNIT(I_TEMP),          & ! [IN]
                          'XY',                       OCEAN_RESTART_OUT_DTYPE,   nohalo=.true.              ) ! [IN]
       call FILEIO_write( OCEAN_SFC_TEMP(:,:),        basename,                  OCEAN_RESTART_OUT_TITLE,   & ! [IN]
                          VAR_NAME(I_SFC_TEMP),       VAR_DESC(I_SFC_TEMP),      VAR_UNIT(I_SFC_TEMP),      & ! [IN]
                          'XY',                       OCEAN_RESTART_OUT_DTYPE,   nohalo=.true.              ) ! [IN]
       call FILEIO_write( OCEAN_SFC_albedo(:,:,I_LW), basename,                  OCEAN_RESTART_OUT_TITLE,   & ! [IN]
                          VAR_NAME(I_SFC_albedo_LW),  VAR_DESC(I_SFC_albedo_LW), VAR_UNIT(I_SFC_albedo_LW), & ! [IN]
                          'XY',                       OCEAN_RESTART_OUT_DTYPE,   nohalo=.true.              ) ! [IN]
       call FILEIO_write( OCEAN_SFC_albedo(:,:,I_SW), basename,                  OCEAN_RESTART_OUT_TITLE,   & ! [IN]
                          VAR_NAME(I_SFC_albedo_SW),  VAR_DESC(I_SFC_albedo_SW), VAR_UNIT(I_SFC_albedo_SW), & ! [IN]
                          'XY',                       OCEAN_RESTART_OUT_DTYPE,   nohalo=.true.              ) ! [IN]
       call FILEIO_write( OCEAN_SFC_Z0M(:,:),         basename,                  OCEAN_RESTART_OUT_TITLE,   & ! [IN]
                          VAR_NAME(I_SFC_Z0M),        VAR_DESC(I_SFC_Z0M),       VAR_UNIT(I_SFC_Z0M),       & ! [IN]
                          'XY',                       OCEAN_RESTART_OUT_DTYPE,   nohalo=.true.              ) ! [IN]
       call FILEIO_write( OCEAN_SFC_Z0H(:,:),         basename,                  OCEAN_RESTART_OUT_TITLE,   & ! [IN]
                          VAR_NAME(I_SFC_Z0H),        VAR_DESC(I_SFC_Z0H),       VAR_UNIT(I_SFC_Z0H),       & ! [IN]
                          'XY',                       OCEAN_RESTART_OUT_DTYPE,   nohalo=.true.              ) ! [IN]
       call FILEIO_write( OCEAN_SFC_Z0E(:,:),         basename,                  OCEAN_RESTART_OUT_TITLE,   & ! [IN]
                          VAR_NAME(I_SFC_Z0E),        VAR_DESC(I_SFC_Z0E),       VAR_UNIT(I_SFC_Z0E),       & ! [IN]
                          'XY',                       OCEAN_RESTART_OUT_DTYPE,   nohalo=.true.              ) ! [IN]
       call FILEIO_write( OCEAN_SFLX_heat(:,:),       basename,                  OCEAN_RESTART_OUT_TITLE,   & ! [IN]
                          VAR_NAME(I_SFLX_heat),      VAR_DESC(I_SFLX_heat),     VAR_UNIT(I_SFLX_heat),     & ! [IN]
                          'XY',                       OCEAN_RESTART_OUT_DTYPE,   nohalo=.true.              ) ! [IN]
       call FILEIO_write( OCEAN_SFLX_prec(:,:),       basename,                  OCEAN_RESTART_OUT_TITLE,   & ! [IN]
                          VAR_NAME(I_SFLX_prec),      VAR_DESC(I_SFLX_prec),     VAR_UNIT(I_SFLX_prec),     & ! [IN]
                          'XY',                       OCEAN_RESTART_OUT_DTYPE,   nohalo=.true.              ) ! [IN]
       call FILEIO_write( OCEAN_SFLX_evap(:,:),       basename,                  OCEAN_RESTART_OUT_TITLE,   & ! [IN]
                          VAR_NAME(I_SFLX_evap),      VAR_DESC(I_SFLX_evap),     VAR_UNIT(I_SFLX_evap),     & ! [IN]
                          'XY',                       OCEAN_RESTART_OUT_DTYPE,   nohalo=.true.              ) ! [IN]

    endif

    return
  end subroutine OCEAN_vars_restart_write

  !-----------------------------------------------------------------------------
  !> History output set for ocean variables
  subroutine OCEAN_vars_history
    use scale_time, only: &
       TIME_DTSEC_OCEAN
    use scale_history, only: &
       HIST_in
    implicit none
    !---------------------------------------------------------------------------

    if ( OCEAN_VARS_CHECKRANGE ) then
       call VALCHECK( OCEAN_TEMP      (IS:IE,JS:JE),      0.0_RP, 1000.0_RP, VAR_NAME(I_TEMP),          &
                     __FILE__, __LINE__ )
       call VALCHECK( OCEAN_SFC_TEMP  (IS:IE,JS:JE),      0.0_RP, 1000.0_RP, VAR_NAME(I_SFC_TEMP),      &
                     __FILE__, __LINE__ )
       call VALCHECK( OCEAN_SFC_albedo(IS:IE,JS:JE,I_LW), 0.0_RP,    2.0_RP, VAR_NAME(I_SFC_albedo_LW), &
                     __FILE__, __LINE__ )
       call VALCHECK( OCEAN_SFC_albedo(IS:IE,JS:JE,I_SW), 0.0_RP,    2.0_RP, VAR_NAME(I_SFC_albedo_SW), &
                     __FILE__, __LINE__ )
       call VALCHECK( OCEAN_SFC_Z0M   (IS:IE,JS:JE),      0.0_RP, 1000.0_RP, VAR_NAME(I_SFC_Z0M),       &
                     __FILE__, __LINE__ )
       call VALCHECK( OCEAN_SFC_Z0H   (IS:IE,JS:JE),      0.0_RP, 1000.0_RP, VAR_NAME(I_SFC_Z0H),       &
                     __FILE__, __LINE__ )
       call VALCHECK( OCEAN_SFC_Z0E   (IS:IE,JS:JE),      0.0_RP, 1000.0_RP, VAR_NAME(I_SFC_Z0E),       &
                     __FILE__, __LINE__ )
    endif

    call HIST_in( OCEAN_TEMP      (:,:),      VAR_NAME(I_TEMP),          VAR_DESC(I_TEMP),          VAR_UNIT(I_TEMP)          )
    call HIST_in( OCEAN_SFC_TEMP  (:,:),      VAR_NAME(I_SFC_TEMP),      VAR_DESC(I_SFC_TEMP),      VAR_UNIT(I_SFC_TEMP)      )
    call HIST_in( OCEAN_SFC_albedo(:,:,I_LW), VAR_NAME(I_SFC_albedo_LW), VAR_DESC(I_SFC_albedo_LW), VAR_UNIT(I_SFC_albedo_LW) )
    call HIST_in( OCEAN_SFC_albedo(:,:,I_SW), VAR_NAME(I_SFC_albedo_SW), VAR_DESC(I_SFC_albedo_SW), VAR_UNIT(I_SFC_albedo_SW) )
    call HIST_in( OCEAN_SFC_Z0M   (:,:),      VAR_NAME(I_SFC_Z0M),       VAR_DESC(I_SFC_Z0M),       VAR_UNIT(I_SFC_Z0M)       )
    call HIST_in( OCEAN_SFC_Z0H   (:,:),      VAR_NAME(I_SFC_Z0H),       VAR_DESC(I_SFC_Z0H),       VAR_UNIT(I_SFC_Z0H)       )
    call HIST_in( OCEAN_SFC_Z0E   (:,:),      VAR_NAME(I_SFC_Z0E),       VAR_DESC(I_SFC_Z0E),       VAR_UNIT(I_SFC_Z0E)       )
    call HIST_in( OCEAN_SFLX_heat (:,:),      VAR_NAME(I_SFLX_heat),     VAR_DESC(I_SFLX_heat),     VAR_UNIT(I_SFLX_heat)     )
    call HIST_in( OCEAN_SFLX_prec (:,:),      VAR_NAME(I_SFLX_prec),     VAR_DESC(I_SFLX_prec),     VAR_UNIT(I_SFLX_prec)     )
    call HIST_in( OCEAN_SFLX_evap (:,:),      VAR_NAME(I_SFLX_evap),     VAR_DESC(I_SFLX_evap),     VAR_UNIT(I_SFLX_evap)     )

    return
  end subroutine OCEAN_vars_history

  !-----------------------------------------------------------------------------
  !> Budget monitor for ocean
  subroutine OCEAN_vars_total
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STAT_total
    implicit none

    real(RP) :: total
    !---------------------------------------------------------------------------

    if ( STATISTICS_checktotal ) then

       call STAT_total( total, OCEAN_TEMP      (:,:),      VAR_NAME(I_TEMP         ) )
       call STAT_total( total, OCEAN_SFC_TEMP  (:,:),      VAR_NAME(I_SFC_TEMP     ) )
       call STAT_total( total, OCEAN_SFC_albedo(:,:,I_LW), VAR_NAME(I_SFC_albedo_LW) )
       call STAT_total( total, OCEAN_SFC_albedo(:,:,I_SW), VAR_NAME(I_SFC_albedo_SW) )
       call STAT_total( total, OCEAN_SFC_Z0M   (:,:),      VAR_NAME(I_SFC_Z0M      ) )
       call STAT_total( total, OCEAN_SFC_Z0H   (:,:),      VAR_NAME(I_SFC_Z0H      ) )
       call STAT_total( total, OCEAN_SFC_Z0E   (:,:),      VAR_NAME(I_SFC_Z0E      ) )

    endif

    return
  end subroutine OCEAN_vars_total

  !-----------------------------------------------------------------------------
  !> Input from External I/O
  subroutine OCEAN_vars_external_in( &
       OCEAN_TEMP_in,       &
       OCEAN_SFC_TEMP_in,   &
       OCEAN_SFC_albedo_in, &
       OCEAN_SFC_Z0M_in,    &
       OCEAN_SFC_Z0H_in,    &
       OCEAN_SFC_Z0E_in     )
    implicit none

    real(RP), intent(in) :: OCEAN_TEMP_in      (IA,JA)
    real(RP), intent(in) :: OCEAN_SFC_TEMP_in  (IA,JA)
    real(RP), intent(in) :: OCEAN_SFC_albedo_in(IA,JA,2)
    real(RP), intent(in) :: OCEAN_SFC_Z0M_in   (IA,JA)
    real(RP), intent(in) :: OCEAN_SFC_Z0H_in   (IA,JA)
    real(RP), intent(in) :: OCEAN_SFC_Z0E_in   (IA,JA)
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** External Input file (ocean) ***'

    OCEAN_TEMP      (:,:)   = OCEAN_TEMP_in      (:,:)
    OCEAN_SFC_TEMP  (:,:)   = OCEAN_SFC_TEMP_in  (:,:)
    OCEAN_SFC_albedo(:,:,:) = OCEAN_SFC_albedo_in(:,:,:)
    OCEAN_SFC_Z0M   (:,:)   = OCEAN_SFC_Z0M_in   (:,:)
    OCEAN_SFC_Z0H   (:,:)   = OCEAN_SFC_Z0H_in   (:,:)
    OCEAN_SFC_Z0E   (:,:)   = OCEAN_SFC_Z0E_in   (:,:)

    call OCEAN_vars_total

    return
  end subroutine OCEAN_vars_external_in

end module mod_ocean_vars
