!-------------------------------------------------------------------------------
!> module LAND Variables
!!
!! @par Description
!!          Container for land variables
!!
!! @author Team SCALE
!! @li      2013-08-31 (T.Yamaura)  [new]
!<
!-------------------------------------------------------------------------------
module mod_land_vars
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_debug
  use scale_grid_index
  use scale_land_grid_index

  use scale_const, only: &
     I_SW  => CONST_I_SW, &
     I_LW  => CONST_I_LW
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
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
  public :: LAND_vars_external_in

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical, public :: LAND_RESTART_OUTPUT = .false. !< output restart file?

  ! prognostic variables
  real(RP), public, allocatable :: LAND_TEMP      (:,:,:) !< temperature of each soil layer [K]
  real(RP), public, allocatable :: LAND_WATER     (:,:,:) !< moisture    of each soil layer [m3/m3]

  ! tendency variables
  real(RP), public, allocatable :: LAND_TEMP_t    (:,:,:) !< tendency of LAND_TEMP
  real(RP), public, allocatable :: LAND_WATER_t   (:,:,:) !< tendency of LAND_WATER

  ! for restart
  real(RP), public, allocatable :: LAND_SFC_TEMP  (:,:)   !< land surface skin temperature [K]
  real(RP), public, allocatable :: LAND_SFC_albedo(:,:,:) !< land surface albedo           [0-1]

  real(RP), public, allocatable :: LAND_PROPERTY  (:,:,:) !< land surface property

  integer,  public, parameter   :: LAND_PROPERTY_nmax = 8
  integer,  public, parameter   :: I_WaterLimit       = 1 ! maximum  soil moisture        [m3/m3]
  integer,  public, parameter   :: I_WaterCritical    = 2 ! critical soil moisture        [m3/m3]
  integer,  public, parameter   :: I_ThermalCond      = 3 ! thermal conductivity for soil [W/m/K]
  integer,  public, parameter   :: I_HeatCapacity     = 4 ! heat capacity        for soil [J/K]
  integer,  public, parameter   :: I_WaterDiff        = 5 ! moisture diffusivity in the soil [m2/s]
  integer,  public, parameter   :: I_Z0M              = 6 ! roughness length for momemtum [m]
  integer,  public, parameter   :: I_Z0H              = 7 ! roughness length for heat     [m]
  integer,  public, parameter   :: I_Z0E              = 8 ! roughness length for moisture [m]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: LAND_param_read

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_LONG),  private :: LAND_RESTART_IN_BASENAME  = ''             !< basename of the restart file
  character(len=H_LONG),  private :: LAND_RESTART_OUT_BASENAME = ''             !< basename of the output file
  character(len=H_MID),   private :: LAND_RESTART_OUT_TITLE    = 'LAND restart' !< title    of the output file
  character(len=H_MID),   private :: LAND_RESTART_OUT_DTYPE    = 'DEFAULT'      !< REAL4 or REAL8

  logical,                private :: LAND_VARS_CHECKRANGE      = .false.

  integer,                private, parameter :: VMAX            = 5 !< number of the variables
  integer,                private, parameter :: I_TEMP          = 1
  integer,                private, parameter :: I_WATER         = 2
  integer,                private, parameter :: I_SFC_TEMP      = 3
  integer,                private, parameter :: I_SFC_albedo_LW = 4
  integer,                private, parameter :: I_SFC_albedo_SW = 5

  character(len=H_SHORT), private            :: VAR_NAME(VMAX) !< name  of the variables
  character(len=H_MID),   private            :: VAR_DESC(VMAX) !< desc. of the variables
  character(len=H_SHORT), private            :: VAR_UNIT(VMAX) !< unit  of the variables

  data VAR_NAME / 'LAND_TEMP',       &
                  'LAND_WATER',      &
                  'LAND_SFC_TEMP',   &
                  'LAND_SFC_ALB_LW', &
                  'LAND_SFC_ALB_SW'  /
  data VAR_DESC / 'temperature at each soil layer',  &
                  'moisture    at each soil layer',  &
                  'land surface skin temperature',   &
                  'land surface albedo (longwave)',  &
                  'land surface albedo (shortwave)'  /
  data VAR_UNIT / 'K',     &
                  'm3/m3', &
                  'K',     &
                  '0-1',   &
                  '0-1'    /

  integer,  private, parameter :: LAND_NUM_IDX = 2 ! # of land indices

  real(RP), private            :: LAND_PROPERTY_table(LAND_NUM_IDX,LAND_PROPERTY_nmax)

  integer,  private              :: LAND_QA_comm
  real(RP), private, allocatable :: work_comm(:,:,:) ! for communication

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine LAND_vars_setup
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    use scale_landuse, only: &
       LANDUSE_index_PFT
    implicit none

    NAMELIST / PARAM_LAND_VARS /  &
       LAND_RESTART_IN_BASENAME,  &
       LAND_RESTART_OUTPUT,       &
       LAND_RESTART_OUT_BASENAME, &
       LAND_RESTART_OUT_TITLE,    &
       LAND_RESTART_OUT_DTYPE,    &
       LAND_VARS_CHECKRANGE

    integer :: ierr
    integer :: i, j, iv, p
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[VARS] / Categ[LAND] / Origin[SCALE-LES]'

    allocate( LAND_TEMP   (LKMAX,IA,JA) )
    allocate( LAND_WATER  (LKMAX,IA,JA) )
    LAND_TEMP   (:,:,:) = UNDEF
    LAND_WATER  (:,:,:) = UNDEF

    allocate( LAND_TEMP_t (LKMAX,IA,JA) )
    allocate( LAND_WATER_t(LKMAX,IA,JA) )
    LAND_TEMP_t (:,:,:) = UNDEF
    LAND_WATER_t(:,:,:) = UNDEF

    allocate( LAND_SFC_TEMP  (IA,JA)   )
    allocate( LAND_SFC_albedo(IA,JA,2) )
    LAND_SFC_TEMP  (:,:)   = UNDEF
    LAND_SFC_albedo(:,:,:) = UNDEF

    LAND_QA_comm = LKMAX &
                 + LKMAX &
                 + 1     &
                 + 2

    allocate( work_comm(IA,JA,LAND_QA_comm) )

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_LAND_VARS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_LAND_VARS. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_LAND_VARS)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** List of prognostic variables (LAND) ***'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A15,A,A32,3(A))') &
               '***       |','VARNAME        ','|', 'DESCRIPTION                     ','[', 'UNIT            ',']'
    do iv = 1, VMAX
       if( IO_L ) write(IO_FID_LOG,'(1x,A,i3,A,A15,A,A32,3(A))') &
                  '*** NO.',iv,'|',VAR_NAME(iv),'|',VAR_DESC(iv),'[',VAR_UNIT(iv),']'
    enddo

    if( IO_L ) write(IO_FID_LOG,*)
    if ( LAND_RESTART_IN_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart input?  : ', trim(LAND_RESTART_IN_BASENAME)
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart input?  : NO'
    endif
    if (       LAND_RESTART_OUTPUT             &
         .AND. LAND_RESTART_OUT_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart output? : ', trim(LAND_RESTART_OUT_BASENAME)
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart output? : NO'
       LAND_RESTART_OUTPUT = .false.
    endif



    call LAND_param_read

    allocate( LAND_PROPERTY(IA,JA,LAND_PROPERTY_nmax) )

    ! tentative, mosaic is off
    do p = 1, LAND_PROPERTY_nmax
    do j = JS, JE
    do i = IS, IE
       LAND_PROPERTY(i,j,p) = LAND_PROPERTY_table( LANDUSE_index_PFT(i,j,1),p )
    enddo
    enddo
    enddo

    do p = 1, LAND_PROPERTY_nmax
       call COMM_vars8( LAND_PROPERTY(:,:,p), p )
    enddo
    do p = 1, LAND_PROPERTY_nmax
       call COMM_wait ( LAND_PROPERTY(:,:,p), p )
    enddo

    return
  end subroutine LAND_vars_setup

  !-----------------------------------------------------------------------------
  !> HALO Communication
  subroutine LAND_vars_fillhalo
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    integer :: k, iq
    !---------------------------------------------------------------------------

    iq = 0
    do k = LKS, LKE
       iq = iq + 1
       work_comm(:,:,iq) = LAND_TEMP (k,:,:)
    enddo
    do k = LKS, LKE
       iq = iq + 1
       work_comm(:,:,iq) = LAND_WATER(k,:,:)
    enddo
    iq = iq + 1
    work_comm(:,:,iq) = LAND_SFC_TEMP  (:,:)
    iq = iq + 1
    work_comm(:,:,iq) = LAND_SFC_albedo(:,:,I_LW)
    iq = iq + 1
    work_comm(:,:,iq) = LAND_SFC_albedo(:,:,I_SW)

    do iq = 1, LAND_QA_comm
      call COMM_vars8( work_comm(:,:,iq), iq )
    enddo
    do iq = 1, LAND_QA_comm
      call COMM_wait ( work_comm(:,:,iq), iq )
    enddo

    iq = 0
    do k = LKS, LKE
       iq = iq + 1
       LAND_TEMP (k,:,:) = work_comm(:,:,iq)
    enddo
    do k = LKS, LKE
       iq = iq + 1
       LAND_WATER(k,:,:) = work_comm(:,:,iq)
    enddo
    iq = iq + 1
    LAND_SFC_TEMP  (:,:)      = work_comm(:,:,iq)
    iq = iq + 1
    LAND_SFC_albedo(:,:,I_LW) = work_comm(:,:,iq)
    iq = iq + 1
    LAND_SFC_albedo(:,:,I_SW) = work_comm(:,:,iq)

    return
  end subroutine LAND_vars_fillhalo

  !-----------------------------------------------------------------------------
  !> Read land restart
  subroutine LAND_vars_restart_read
    use scale_fileio, only: &
       FILEIO_read
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input restart file (LAND) ***'

    if ( LAND_RESTART_IN_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(LAND_RESTART_IN_BASENAME)

       call FILEIO_read( LAND_TEMP (:,:,:),                                                & ! [OUT]
                         LAND_RESTART_IN_BASENAME, VAR_NAME(I_TEMP),        'Land', step=1 ) ! [IN]
       call FILEIO_read( LAND_WATER(:,:,:),                                                & ! [OUT]
                         LAND_RESTART_IN_BASENAME, VAR_NAME(I_WATER),       'Land', step=1 ) ! [IN]
       call FILEIO_read( LAND_SFC_TEMP(:,:),                                               & ! [OUT]
                         LAND_RESTART_IN_BASENAME, VAR_NAME(I_SFC_TEMP),      'XY', step=1 ) ! [IN]
       call FILEIO_read( LAND_SFC_albedo(:,:,I_LW),                                        & ! [OUT]
                         LAND_RESTART_IN_BASENAME, VAR_NAME(I_SFC_albedo_LW), 'XY', step=1 ) ! [IN]
       call FILEIO_read( LAND_SFC_albedo(:,:,I_SW),                                        & ! [OUT]
                         LAND_RESTART_IN_BASENAME, VAR_NAME(I_SFC_albedo_SW), 'XY', step=1 ) ! [IN]

       call LAND_vars_fillhalo

       call LAND_vars_total
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** restart file for land is not specified.'
    endif

    return
  end subroutine LAND_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Write land restart
  subroutine LAND_vars_restart_write
    use scale_time, only: &
       TIME_gettimelabel
    use scale_fileio, only: &
       FILEIO_write
    use mod_cpl_admin, only: &
       CPL_sw_AtmLnd
    use mod_cpl_vars, only: &
       CPL_getLnd_restart
    implicit none

    character(len=15)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    if ( LAND_RESTART_OUT_BASENAME /= '' ) then

       if ( CPL_sw_AtmLnd ) then
          call CPL_getLnd_restart( LAND_SFC_TEMP  (:,:),  & ! [OUT]
                                   LAND_SFC_albedo(:,:,:) ) ! [OUT]
       endif

       call TIME_gettimelabel( timelabel )
       write(basename,'(A,A,A)') trim(LAND_RESTART_OUT_BASENAME), '_', trim(timelabel)

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Output restart file (LAND) ***'
       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(basename)

       call LAND_vars_fillhalo

       call LAND_vars_total

       call FILEIO_write( LAND_TEMP    (:,:,:),      basename,                  LAND_RESTART_OUT_TITLE,    & ! [IN]
                          VAR_NAME(I_TEMP),          VAR_DESC(I_TEMP),          VAR_UNIT(I_TEMP),          & ! [IN]
                          'Land',                                               LAND_RESTART_OUT_DTYPE     ) ! [IN]
       call FILEIO_write( LAND_WATER   (:,:,:),      basename,                  LAND_RESTART_OUT_TITLE,    & ! [IN]
                          VAR_NAME(I_WATER),         VAR_DESC(I_WATER),         VAR_UNIT(I_WATER),         & ! [IN]
                          'Land',                                               LAND_RESTART_OUT_DTYPE     ) ! [IN]
       call FILEIO_write( LAND_SFC_TEMP(:,:),        basename,                  LAND_RESTART_OUT_TITLE,    & ! [IN]
                          VAR_NAME(I_SFC_TEMP),      VAR_DESC(I_SFC_TEMP),      VAR_UNIT(I_SFC_TEMP),      & ! [IN]
                          'XY',                                                 LAND_RESTART_OUT_DTYPE     ) ! [IN]
       call FILEIO_write( LAND_SFC_albedo(:,:,I_LW), basename,                  LAND_RESTART_OUT_TITLE,    & ! [IN]
                          VAR_NAME(I_SFC_albedo_LW), VAR_DESC(I_SFC_albedo_LW), VAR_UNIT(I_SFC_albedo_LW), & ! [IN]
                          'XY',                                                 LAND_RESTART_OUT_DTYPE     ) ! [IN]
       call FILEIO_write( LAND_SFC_albedo(:,:,I_SW), basename,                  LAND_RESTART_OUT_TITLE,    & ! [IN]
                          VAR_NAME(I_SFC_albedo_SW), VAR_DESC(I_SFC_albedo_SW), VAR_UNIT(I_SFC_albedo_SW), & ! [IN]
                          'XY',                                                 LAND_RESTART_OUT_DTYPE     ) ! [IN]

    endif

    return
  end subroutine LAND_vars_restart_write

  !-----------------------------------------------------------------------------
  !> History output set for land variables
  subroutine LAND_vars_history
    use scale_time, only: &
       TIME_DTSEC_LAND
    use scale_history, only: &
       HIST_in
    implicit none
    !---------------------------------------------------------------------------

    if ( LAND_VARS_CHECKRANGE ) then
       call VALCHECK( LAND_TEMP (:,:,:), 0.0_RP, 1000.0_RP, VAR_NAME(I_TEMP) , __FILE__, __LINE__ )
       call VALCHECK( LAND_WATER(:,:,:), 0.0_RP, 1000.0_RP, VAR_NAME(I_WATER), __FILE__, __LINE__ )
    endif

    call HIST_in( LAND_TEMP (:,:,:), VAR_NAME(I_TEMP),  VAR_DESC(I_TEMP),  VAR_UNIT(I_TEMP),  TIME_DTSEC_LAND, zdim='land' )
    call HIST_in( LAND_WATER(:,:,:), VAR_NAME(I_WATER), VAR_DESC(I_WATER), VAR_UNIT(I_WATER), TIME_DTSEC_LAND, zdim='land' )

    return
  end subroutine LAND_vars_history

  !-----------------------------------------------------------------------------
  !> Budget monitor for land
  subroutine LAND_vars_total
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STAT_total
    implicit none

    real(RP) :: total

    character(len=2) :: sk
    integer          :: k
    !---------------------------------------------------------------------------

    if ( STATISTICS_checktotal ) then

       do k = LKS, LKE
          write(sk,'(I2.2)') k

          call STAT_total( total, LAND_TEMP (k,:,:), trim(VAR_NAME(I_TEMP) )//sk )
          call STAT_total( total, LAND_WATER(k,:,:), trim(VAR_NAME(I_WATER))//sk )
       enddo

       call STAT_total( total, LAND_SFC_TEMP  (:,:),      VAR_NAME(I_SFC_TEMP     ) )
       call STAT_total( total, LAND_SFC_albedo(:,:,I_LW), VAR_NAME(I_SFC_albedo_LW) )
       call STAT_total( total, LAND_SFC_albedo(:,:,I_SW), VAR_NAME(I_SFC_albedo_SW) )

    endif

    return
  end subroutine LAND_vars_total

  !-----------------------------------------------------------------------------
  !> Input from External I/O
  subroutine LAND_vars_external_in( &
      LAND_TEMP_in,      &
      LAND_WATER_in,     &
      LAND_SFC_TEMP_in,  &
      LAND_SFC_albedo_in )
    implicit none

    real(RP), intent(in) :: LAND_TEMP_in (:,:,:)
    real(RP), intent(in) :: LAND_WATER_in(:,:,:)
    real(RP), intent(in) :: LAND_SFC_TEMP_in  (IA,JA)
    real(RP), intent(in) :: LAND_SFC_albedo_in(IA,JA,2)
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** External Input (land) ***'

    LAND_TEMP      (:,:,:) = LAND_TEMP_in      (:,:,:)
    LAND_WATER     (:,:,:) = LAND_WATER_in     (:,:,:)
    LAND_SFC_TEMP  (:,:)   = LAND_SFC_TEMP_in  (:,:)
    LAND_SFC_albedo(:,:,:) = LAND_SFC_albedo_in(:,:,:)

    call LAND_vars_fillhalo

    call LAND_vars_total

    return
  end subroutine LAND_vars_external_in

  !-----------------------------------------------------------------------------
  !> Budget monitor for land
  subroutine LAND_param_read
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       CONST_UNDEF
    implicit none

    integer                :: index
    character(len=H_SHORT) :: description
    real(RP)               :: STRGMAX
    real(RP)               :: STRGCRT
    real(RP)               :: TCS
    real(RP)               :: HCS
    real(RP)               :: DFW
    real(RP)               :: Z0M
    real(RP)               :: Z0H
    real(RP)               :: Z0E

    NAMELIST / PARAM_LAND_DATA / &
       index,       &
       description, &
       STRGMAX,     &
       STRGCRT,     &
       TCS,         &
       HCS,         &
       DFW,         &
       Z0M,         &
       Z0H,         &
       Z0E

    integer :: n
    integer :: ierr
    !---------------------------------------------------------------------------

    LAND_PROPERTY_table(:,:) = CONST_UNDEF

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** [LAND ] vegetation parameters'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,11(1x,A))') '***        ',  &
                                                   ' description', &
                                                   ' Max Stg.', &
                                                   ' CRT Stg.', &
                                                   ' T condu.', &
                                                   ' H capac.', &
                                                   ' DFC Wat.', &
                                                   '    Z0(m)', &
                                                   '    Z0(h)', &
                                                   '    Z0(e)'

    !--- read namelist
    rewind(IO_FID_CONF)
    do n = 1, LAND_NUM_IDX
       ! undefined roughness length
       Z0H = -1.0_RP
       Z0E = -1.0_RP

       read(IO_FID_CONF,nml=PARAM_LAND_DATA,iostat=ierr)
       if ( ierr < 0 ) then !--- no more data
          exit
       elseif( ierr > 0 ) then !--- fatal error
          write(*,*) 'xxx Not appropriate names in namelist PARAM_LAND. Check!'
          call PRC_MPIstop
       endif

       if( Z0H < 0.0_RP ) then
         Z0H = Z0M / 7.4_RP ! defined by Garratt and Francey (1978)
       endif
       if( Z0E < 0.0_RP ) then
         Z0E = Z0M / 7.4_RP ! defined by Garratt and Francey (1978)
       endif

       LAND_PROPERTY_table(index,I_WaterLimit   ) = STRGMAX
       LAND_PROPERTY_table(index,I_WaterCritical) = STRGCRT
       LAND_PROPERTY_table(index,I_ThermalCond  ) = TCS
       LAND_PROPERTY_table(index,I_HeatCapacity ) = HCS
       LAND_PROPERTY_table(index,I_WaterDiff    ) = DFW
       LAND_PROPERTY_table(index,I_Z0M          ) = Z0M
       LAND_PROPERTY_table(index,I_Z0H          ) = Z0H
       LAND_PROPERTY_table(index,I_Z0E          ) = Z0E

       if( IO_L ) write(IO_FID_LOG,'(1x,A8,I3,1x,A12,8(1x,F9.2))') &
                                     '*** IDX=', index, &
                                     trim(description), &
                                     STRGMAX, &
                                     STRGCRT, &
                                     TCS,     &
                                     HCS,     &
                                     DFW,     &
                                     Z0M,     &
                                     Z0H,     &
                                     Z0E
    enddo

    return
  end subroutine LAND_param_read

end module mod_land_vars
