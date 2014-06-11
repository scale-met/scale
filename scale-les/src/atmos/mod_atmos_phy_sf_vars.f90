!-------------------------------------------------------------------------------
!> module ATMOSPHERIC Surface Variables
!!
!! @par Description
!!          Container for mod_atmos_phy_sf
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2012-03-27 (H.Yashiro)  [new]
!! @li      2013-08-31 (T.Yamaura)  [mod]
!<
!-------------------------------------------------------------------------------
module mod_atmos_phy_sf_vars
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_SF_vars_setup
  public :: ATMOS_PHY_SF_vars_fillhalo
  public :: ATMOS_PHY_SF_vars_restart_read
  public :: ATMOS_PHY_SF_vars_restart_write

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical,  public :: ATMOS_PHY_SF_RESTART_OUTPUT = .false. !< output restart file?

  real(RP), public, allocatable :: ATMOS_PHY_SF_DENS_t(:,:)   ! tendency DENS [    kg/m3/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_MOMZ_t(:,:)   ! tendency MOMZ [m/s*kg/m3/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_MOMX_t(:,:)   ! tendency MOMX [m/s*kg/m3/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_MOMY_t(:,:)   ! tendency MOMY [m/s*kg/m3/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_RHOT_t(:,:)   ! tendency RHOT [K  *kg/m3/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_QTRC_t(:,:,:) ! tendency QTRC [    kg/kg/s]

  ! surface flux (upward positive)
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_MW   (:,:)   ! z-momentum flux (area center) [m/s*kg/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_MU   (:,:)   ! x-momentum flux (area center) [m/s*kg/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_MV   (:,:)   ! y-momentum flux (area center) [m/s*kg/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_SH   (:,:)   ! sensible heat flux [J/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_LH   (:,:)   ! latent   heat flux [J/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_QTRC (:,:,:) ! tracer mass flux [kg/m2/s]

  real(RP), public, allocatable :: ATMOS_PHY_SF_SFC_DENS  (:,:)   ! surface atmosphere density  [kg/m3]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFC_PRES  (:,:)   ! surface atmosphere pressure [Pa]

  real(RP), public, allocatable :: ATMOS_PHY_SF_SFC_TEMP  (:,:)   ! surface skin temperature             [K]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFC_albedo(:,:,:) ! surface albedo, land only            [0-1]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFC_Z0    (:,:)   ! surface roughness length, ocean only [m]

!  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_QEMIS(:,:,:) ! tracer emission   flux [kg/m2/s]
!  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_QDEP (:,:,:) ! tracer deposition flux [kg/m2/s]
!  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_VDEP (:,:,:) ! tracer deposition velocity [m/s]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_LONG),  private :: ATMOS_PHY_SF_RESTART_IN_BASENAME  = ''                     !< basename of the restart file
  character(len=H_LONG),  private :: ATMOS_PHY_SF_RESTART_OUT_BASENAME = ''                     !< basename of the output file
  character(len=H_MID),   private :: ATMOS_PHY_SF_RESTART_OUT_TITLE    = 'ATMOS_PHY_SF restart' !< title    of the output file
  character(len=H_MID),   private :: ATMOS_PHY_SF_RESTART_OUT_DTYPE    = 'DEFAULT'              !< REAL4 or REAL8

  integer,                private, parameter :: VMAX = 4       !< number of the variables
  integer,                private, parameter :: I_SFC_TEMP      = 1
  integer,                private, parameter :: I_SFC_albedo_LW = 2
  integer,                private, parameter :: I_SFC_albedo_SW = 3
  integer,                private, parameter :: I_SFC_Z0        = 4

  character(len=H_SHORT), private            :: VAR_NAME(VMAX) !< name  of the variables
  character(len=H_MID),   private            :: VAR_DESC(VMAX) !< desc. of the variables
  character(len=H_SHORT), private            :: VAR_UNIT(VMAX) !< unit  of the variables

  data VAR_NAME / 'SFC_TEMP', &
                  'ALB_LW',   &
                  'ALB_SW',   &
                  'SFC_Z0'    /

  data VAR_DESC / 'surface skin temperature',     &
                  'surface albedo for longwave',  &
                  'surface albedo for shortwave', &
                  'surface roughness length'      /

  data VAR_UNIT / 'K',   &
                  '0-1', &
                  '0-1', &
                  'm'    /

  real(RP), private :: ATMOS_PHY_SF_DEFAULT_SFC_TEMP
  real(RP), private :: ATMOS_PHY_SF_DEFAULT_SFC_albedo

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_SF_vars_setup
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    implicit none

    NAMELIST / PARAM_ATMOS_PHY_SF_VARS / &
       ATMOS_PHY_SF_RESTART_IN_BASENAME,  &
       ATMOS_PHY_SF_RESTART_OUTPUT,       &
       ATMOS_PHY_SF_RESTART_OUT_BASENAME, &
       ATMOS_PHY_SF_RESTART_OUT_TITLE,    &
       ATMOS_PHY_SF_RESTART_OUT_DTYPE,    &
       ATMOS_PHY_SF_DEFAULT_SFC_TEMP,     &
       ATMOS_PHY_SF_DEFAULT_SFC_albedo

    integer :: ierr
    integer :: iv
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[VARS] / Categ[ATMOS_PHY_SF] / Origin[SCALE-LES]'

    allocate( ATMOS_PHY_SF_DENS_t(IA,JA)    )
    allocate( ATMOS_PHY_SF_MOMZ_t(IA,JA)    )
    allocate( ATMOS_PHY_SF_MOMX_t(IA,JA)    )
    allocate( ATMOS_PHY_SF_MOMY_t(IA,JA)    )
    allocate( ATMOS_PHY_SF_RHOT_t(IA,JA)    )
    allocate( ATMOS_PHY_SF_QTRC_t(IA,JA,QA) )
    ATMOS_PHY_SF_DENS_t(:,:)   = UNDEF
    ATMOS_PHY_SF_MOMZ_t(:,:)   = UNDEF
    ATMOS_PHY_SF_MOMX_t(:,:)   = UNDEF
    ATMOS_PHY_SF_MOMY_t(:,:)   = UNDEF
    ATMOS_PHY_SF_RHOT_t(:,:)   = UNDEF
    ATMOS_PHY_SF_QTRC_t(:,:,:) = UNDEF

    allocate( ATMOS_PHY_SF_SFLX_MW   (IA,JA)    )
    allocate( ATMOS_PHY_SF_SFLX_MU   (IA,JA)    )
    allocate( ATMOS_PHY_SF_SFLX_MV   (IA,JA)    )
    allocate( ATMOS_PHY_SF_SFLX_SH   (IA,JA)    )
    allocate( ATMOS_PHY_SF_SFLX_LH   (IA,JA)    )
    allocate( ATMOS_PHY_SF_SFLX_QTRC (IA,JA,QA) )
    ATMOS_PHY_SF_SFLX_MW  (:,:)   = UNDEF
    ATMOS_PHY_SF_SFLX_MU  (:,:)   = UNDEF
    ATMOS_PHY_SF_SFLX_MV  (:,:)   = UNDEF
    ATMOS_PHY_SF_SFLX_SH  (:,:)   = UNDEF
    ATMOS_PHY_SF_SFLX_LH  (:,:)   = UNDEF
    ATMOS_PHY_SF_SFLX_QTRC(:,:,:) = UNDEF

    allocate( ATMOS_PHY_SF_SFC_DENS  (IA,JA)    )
    allocate( ATMOS_PHY_SF_SFC_PRES  (IA,JA)    )
    allocate( ATMOS_PHY_SF_SFC_TEMP  (IA,JA)    )
    allocate( ATMOS_PHY_SF_SFC_albedo(IA,JA,2)  )
    allocate( ATMOS_PHY_SF_SFC_Z0    (IA,JA)    )
    ATMOS_PHY_SF_SFC_DENS  (:,:)   = UNDEF
    ATMOS_PHY_SF_SFC_PRES  (:,:)   = UNDEF
    ATMOS_PHY_SF_SFC_TEMP  (:,:)   = UNDEF
    ATMOS_PHY_SF_SFC_albedo(:,:,:) = UNDEF
    ATMOS_PHY_SF_SFC_Z0    (:,:)   = UNDEF

    ATMOS_PHY_SF_DEFAULT_SFC_TEMP   = UNDEF
    ATMOS_PHY_SF_DEFAULT_SFC_albedo = UNDEF

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_SF_VARS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_SF_VARS. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_SF_VARS)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** [ATMOS_PHY_SF] prognostic/diagnostic variables'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A15,A,A32,3(A))') &
               '***       |','VARNAME        ','|', 'DESCRIPTION                     ','[', 'UNIT            ',']'
    do iv = 1, VMAX
       if( IO_L ) write(IO_FID_LOG,'(1x,A,i3,A,A15,A,A32,3(A))') &
                  '*** NO.',iv,'|',VAR_NAME(iv),'|',VAR_DESC(iv),'[',VAR_UNIT(iv),']'
    enddo

    if( IO_L ) write(IO_FID_LOG,*)
    if ( ATMOS_PHY_SF_RESTART_IN_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart input?  : ', trim(ATMOS_PHY_SF_RESTART_IN_BASENAME)
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart input?  : NO'
    endif
    if (       ATMOS_PHY_SF_RESTART_OUTPUT             &
         .AND. ATMOS_PHY_SF_RESTART_OUT_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart output? : ', trim(ATMOS_PHY_SF_RESTART_OUT_BASENAME)
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart output? : NO'
       ATMOS_PHY_SF_RESTART_OUTPUT = .false.
    endif

    return
  end subroutine ATMOS_PHY_SF_vars_setup

  !-----------------------------------------------------------------------------
  !> HALO Communication
  subroutine ATMOS_PHY_SF_vars_fillhalo
    use scale_const, only: &
       I_SW => CONST_I_SW, &
       I_LW => CONST_I_LW
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none
    !---------------------------------------------------------------------------

    call COMM_vars8( ATMOS_PHY_SF_SFC_TEMP  (:,:),      1 )
    call COMM_vars8( ATMOS_PHY_SF_SFC_albedo(:,:,I_LW), 2 )
    call COMM_vars8( ATMOS_PHY_SF_SFC_albedo(:,:,I_SW), 3 )
    call COMM_vars8( ATMOS_PHY_SF_SFC_Z0    (:,:),      4 )
    call COMM_wait ( ATMOS_PHY_SF_SFC_TEMP  (:,:),      1 )
    call COMM_wait ( ATMOS_PHY_SF_SFC_albedo(:,:,I_LW), 2 )
    call COMM_wait ( ATMOS_PHY_SF_SFC_albedo(:,:,I_SW), 3 )
    call COMM_wait ( ATMOS_PHY_SF_SFC_Z0    (:,:),      4 )

    return
  end subroutine ATMOS_PHY_SF_vars_fillhalo

  !-----------------------------------------------------------------------------
  !> Read restart
  subroutine ATMOS_PHY_SF_vars_restart_read
    use scale_const, only: &
       I_SW => CONST_I_SW, &
       I_LW => CONST_I_LW
    use scale_fileio, only: &
       FILEIO_read
    use scale_statistics, only: &
       STAT_total
    implicit none

    character(len=H_SHORT) :: TRC_NAME

    real(RP) :: total
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input restart file (ATMOS_PHY_SF) ***'

    if ( ATMOS_PHY_SF_RESTART_IN_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(ATMOS_PHY_SF_RESTART_IN_BASENAME)

       call FILEIO_read( ATMOS_PHY_SF_SFC_TEMP  (:,:),                               & ! [OUT]
                         ATMOS_PHY_SF_RESTART_IN_BASENAME, VAR_NAME(1), 'XY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_SF_SFC_albedo(:,:,I_LW),                          & ! [OUT]
                         ATMOS_PHY_SF_RESTART_IN_BASENAME, VAR_NAME(2), 'XY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_SF_SFC_albedo(:,:,I_SW),                          & ! [OUT]
                         ATMOS_PHY_SF_RESTART_IN_BASENAME, VAR_NAME(3), 'XY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_SF_SFC_Z0    (:,:),                               & ! [OUT]
                         ATMOS_PHY_SF_RESTART_IN_BASENAME, VAR_NAME(4), 'XY', step=1 ) ! [IN]

       call ATMOS_PHY_SF_vars_fillhalo

       call STAT_total( total, ATMOS_PHY_SF_SFC_TEMP  (:,:),      VAR_NAME(1) )
       call STAT_total( total, ATMOS_PHY_SF_SFC_albedo(:,:,I_LW), VAR_NAME(2) )
       call STAT_total( total, ATMOS_PHY_SF_SFC_albedo(:,:,I_SW), VAR_NAME(3) )
       call STAT_total( total, ATMOS_PHY_SF_SFC_Z0    (:,:),      VAR_NAME(4) )
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** restart file for ATMOS_PHY_SF is not specified.'
       ATMOS_PHY_SF_SFC_TEMP  (:,:)   = ATMOS_PHY_SF_DEFAULT_SFC_TEMP
       ATMOS_PHY_SF_SFC_albedo(:,:,:) = ATMOS_PHY_SF_DEFAULT_SFC_albedo
    endif

    return
  end subroutine ATMOS_PHY_SF_vars_restart_read


  !-----------------------------------------------------------------------------
  !> Write restart
  subroutine ATMOS_PHY_SF_vars_restart_write
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
    !---------------------------------------------------------------------------

    if ( ATMOS_PHY_SF_RESTART_OUT_BASENAME /= '' ) then

       call TIME_gettimelabel( timelabel )
       write(basename,'(A,A,A)') trim(ATMOS_PHY_SF_RESTART_OUT_BASENAME), '_', trim(timelabel)

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Output restart file (ATMOS_PHY_SF) ***'
       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(basename)

       call FILEIO_write( ATMOS_PHY_SF_SFC_TEMP  (:,:),      basename, ATMOS_PHY_SF_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(1), VAR_DESC(1), VAR_UNIT(1), 'XY', ATMOS_PHY_SF_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( ATMOS_PHY_SF_SFC_albedo(:,:,I_LW), basename, ATMOS_PHY_SF_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(2), VAR_DESC(2), VAR_UNIT(2), 'XY', ATMOS_PHY_SF_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( ATMOS_PHY_SF_SFC_albedo(:,:,I_SW), basename, ATMOS_PHY_SF_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(3), VAR_DESC(3), VAR_UNIT(3), 'XY', ATMOS_PHY_SF_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( ATMOS_PHY_SF_SFC_Z0    (:,:),      basename, ATMOS_PHY_SF_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(4), VAR_DESC(4), VAR_UNIT(4), 'XY', ATMOS_PHY_SF_RESTART_OUT_DTYPE  ) ! [IN]
    endif

    return
  end subroutine ATMOS_PHY_SF_vars_restart_write

end module mod_atmos_phy_sf_vars
