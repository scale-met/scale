!-------------------------------------------------------------------------------
!> module URBAN Variables
!!
!! @par Description
!!          Container for urban variables
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
module mod_urban_vars
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_debug
  use scale_grid_index
  use scale_urban_grid_index
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
  public :: URBAN_vars_setup
  public :: URBAN_vars_fillhalo
  public :: URBAN_vars_restart_read
  public :: URBAN_vars_restart_write
  public :: URBAN_vars_history
  public :: URBAN_vars_total

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical, public :: URBAN_sw_restart

  ! prognostic variables
  real(RP), public, allocatable :: TR_URB (:,:)   ! Surface temperature of roof [K]
  real(RP), public, allocatable :: TB_URB (:,:)   ! Surface temperature of building [K]
  real(RP), public, allocatable :: TG_URB (:,:)   ! Surface temperature of ground [K]
  real(RP), public, allocatable :: TC_URB (:,:)   ! Diagnostic canopy air temperature [K]
  real(RP), public, allocatable :: QC_URB (:,:)   ! Diagnostic canopy humidity [-]
  real(RP), public, allocatable :: UC_URB (:,:)   ! Diagnostic canopy wind [m/s]
  real(RP), public, allocatable :: TRL_URB(:,:,:) ! temperature in layer of roof [K]
  real(RP), public, allocatable :: TBL_URB(:,:,:) ! temperature in layer of building [K]
  real(RP), public, allocatable :: TGL_URB(:,:,:) ! temperature in layer of ground [K]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical,                private :: URBAN_RESTART_OUTPUT       = .false.         !< output restart file?
  character(len=H_LONG),  private :: URBAN_RESTART_IN_BASENAME  = ''              !< basename of the restart file
  character(len=H_LONG),  private :: URBAN_RESTART_OUT_BASENAME = ''              !< basename of the output file
  character(len=H_MID),   private :: URBAN_RESTART_OUT_TITLE    = 'URBAN restart' !< title    of the output file
  character(len=H_MID),   private :: URBAN_RESTART_OUT_DTYPE    = 'DEFAULT'       !< REAL4 or REAL8

  logical,                private :: URBAN_VARS_CHECKRANGE      = .false.

  integer,                private, parameter :: VMAX = 9
  integer,                private, parameter :: I_TR_URB  = 1
  integer,                private, parameter :: I_TB_URB  = 2
  integer,                private, parameter :: I_TG_URB  = 3
  integer,                private, parameter :: I_TC_URB  = 4
  integer,                private, parameter :: I_QC_URB  = 5
  integer,                private, parameter :: I_UC_URB  = 6
  integer,                private, parameter :: I_TRL_URB = 7
  integer,                private, parameter :: I_TBL_URB = 8
  integer,                private, parameter :: I_TGL_URB = 9

  character(len=H_SHORT), private            :: VAR_NAME(VMAX) !< name  of the urban variables
  character(len=H_MID),   private            :: VAR_DESC(VMAX) !< desc. of the urban variables
  character(len=H_SHORT), private            :: VAR_UNIT(VMAX) !< unit  of the urban variables

  data VAR_NAME / 'TR_URB' , &
                  'TB_URB' , &
                  'TG_URB' , &
                  'TC_URB' , &
                  'QC_URB' , &
                  'UC_URB' , &
                  'TRL_URB', &
                  'TBL_URB', &
                  'TGL_URB'  /

  data VAR_DESC / 'Surface temperature of roof',       &
                  'Surface temperature of building',   &
                  'Surface temperature of ground',     &
                  'Diagnostic canopy air temperature', &
                  'Diagnostic canopy humidity',        &
                  'Diagnostic canopy wind',            &
                  'temperature in layer of roof',      &
                  'temperature in layer of building',  &
                  'temperature in layer of ground'     /

  data VAR_UNIT / 'K',   &
                  'K',   &
                  'K',   &
                  'K',   &
                  '-',   &
                  'm/s', &
                  'K',   &
                  'K',   &
                  'K'    /

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine URBAN_vars_setup
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    implicit none

    NAMELIST / PARAM_URBAN_VARS /  &
       URBAN_RESTART_IN_BASENAME,  &
       URBAN_RESTART_OUTPUT,       &
       URBAN_RESTART_OUT_BASENAME, &
       URBAN_RESTART_OUT_TITLE,    &
       URBAN_RESTART_OUT_DTYPE,    &
       URBAN_VARS_CHECKRANGE

    integer :: ierr
    integer :: ip
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[VARS] / Categ[URBAN] / Origin[SCALE-LES]'

    ! allocate arrays
    allocate( TR_URB(IA,JA) )
    allocate( TB_URB(IA,JA) )
    allocate( TG_URB(IA,JA) )
    allocate( TC_URB(IA,JA) )
    allocate( QC_URB(IA,JA) )
    allocate( UC_URB(IA,JA) )
    TR_URB (:,:) = UNDEF
    TB_URB (:,:) = UNDEF
    TG_URB (:,:) = UNDEF
    TC_URB (:,:) = UNDEF
    QC_URB (:,:) = UNDEF
    UC_URB (:,:) = UNDEF

    allocate( TRL_URB(UKS:UKE,IA,JA) )
    allocate( TBL_URB(UKS:UKE,IA,JA) )
    allocate( TGL_URB(UKS:UKE,IA,JA) )
    TRL_URB(:,:,:) = UNDEF
    TBL_URB(:,:,:) = UNDEF
    TGL_URB(:,:,:) = UNDEF

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_URBAN_VARS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_URBAN_VARS. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_URBAN_VARS)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** List of prognostic variables (URBAN) ***'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A16,A,A32,3(A))') &
               '***       |','         VARNAME','|', 'DESCRIPTION                     ','[', 'UNIT            ',']'
    do ip = 1, VMAX
       if( IO_L ) write(IO_FID_LOG,'(1x,A,i3,A,A16,A,A32,3(A))') &
                  '*** NO.',ip,'|',trim(VAR_NAME(ip)),'|', VAR_DESC(ip),'[', VAR_UNIT(ip),']'
    enddo

    if( IO_L ) write(IO_FID_LOG,*)
    if ( URBAN_RESTART_IN_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart input?  : ', trim(URBAN_RESTART_IN_BASENAME)
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart input?  : NO'
    endif
    if (       URBAN_RESTART_OUTPUT             &
         .AND. URBAN_RESTART_OUT_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart output? : ', trim(URBAN_RESTART_OUT_BASENAME)
       URBAN_sw_restart = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart output? : NO'
       URBAN_RESTART_OUTPUT = .false.
       URBAN_sw_restart = .false.
    endif

    return
  end subroutine URBAN_vars_setup

  !-----------------------------------------------------------------------------
  !> HALO Communication
  subroutine URBAN_vars_fillhalo
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    integer :: k

    real(RP) :: tmp1(IA,JA,UKS:UKE)
    real(RP) :: tmp2(IA,JA,UKS:UKE)
    real(RP) :: tmp3(IA,JA,UKS:UKE)
    !---------------------------------------------------------------------------

    do k = UKS, UKE
      tmp1(:,:,k) = TRL_URB(k,:,:)
      tmp2(:,:,k) = TBL_URB(k,:,:)
      tmp3(:,:,k) = TGL_URB(k,:,:)
    end do

    do k = UKS, UKE
      call COMM_vars8( tmp1(:,:,k), k       )
      call COMM_vars8( tmp2(:,:,k), k+UKE   )
      call COMM_vars8( tmp3(:,:,k), k+UKE*2 )
    end do

    call COMM_vars8( TR_URB(:,:), 1+UKE*3 )
    call COMM_vars8( TB_URB(:,:), 2+UKE*3 )
    call COMM_vars8( TG_URB(:,:), 3+UKE*3 )
    call COMM_vars8( TC_URB(:,:), 4+UKE*3 )
    call COMM_vars8( QC_URB(:,:), 5+UKE*3 )
    call COMM_vars8( UC_URB(:,:), 6+UKE*3 )

    do k = UKS, UKE
      call COMM_wait ( tmp1(:,:,k), k       )
      call COMM_wait ( tmp2(:,:,k), k+UKE   )
      call COMM_wait ( tmp3(:,:,k), k+UKE*2 )
    end do

    call COMM_wait ( TR_URB(:,:), 1+UKE*3 )
    call COMM_wait ( TB_URB(:,:), 2+UKE*3 )
    call COMM_wait ( TG_URB(:,:), 3+UKE*3 )
    call COMM_wait ( TC_URB(:,:), 4+UKE*3 )
    call COMM_wait ( QC_URB(:,:), 5+UKE*3 )
    call COMM_wait ( UC_URB(:,:), 6+UKE*3 )

    do k = UKS, UKE
      TRL_URB(k,:,:) = tmp1(:,:,k)
      TBL_URB(k,:,:) = tmp2(:,:,k)
      TGL_URB(k,:,:) = tmp3(:,:,k)
    end do

    return
  end subroutine URBAN_vars_fillhalo

  !-----------------------------------------------------------------------------
  !> Read urban restart
  subroutine URBAN_vars_restart_read
    use scale_fileio, only: &
       FILEIO_read
    implicit none

    integer :: i, j, v
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input restart file (URBAN) ***'

    if ( URBAN_RESTART_IN_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(URBAN_RESTART_IN_BASENAME)

       call FILEIO_read( TR_URB(:,:),                                      & ! [OUT]
                         URBAN_RESTART_IN_BASENAME, 'TR_URB', 'XY', step=1 ) ! [IN]
       call FILEIO_read( TB_URB(:,:),                                      & ! [OUT]
                         URBAN_RESTART_IN_BASENAME, 'TB_URB', 'XY', step=1 ) ! [IN]
       call FILEIO_read( TG_URB(:,:),                                      & ! [OUT]
                         URBAN_RESTART_IN_BASENAME, 'TG_URB', 'XY', step=1 ) ! [IN]
       call FILEIO_read( TC_URB(:,:),                                      & ! [OUT]
                         URBAN_RESTART_IN_BASENAME, 'TC_URB', 'XY', step=1 ) ! [IN]
       call FILEIO_read( QC_URB(:,:),                                      & ! [OUT]
                         URBAN_RESTART_IN_BASENAME, 'QC_URB', 'XY', step=1 ) ! [IN]
       call FILEIO_read( UC_URB(:,:),                                      & ! [OUT]
                         URBAN_RESTART_IN_BASENAME, 'UC_URB', 'XY', step=1 ) ! [IN]

       call FILEIO_read( TRL_URB(:,:,:),                                       & ! [OUT]
                         URBAN_RESTART_IN_BASENAME, 'TRL_URB', 'Urban', step=1 ) ! [IN]
       call FILEIO_read( TBL_URB(:,:,:),                                       & ! [OUT]
                         URBAN_RESTART_IN_BASENAME, 'TBL_URB', 'Urban', step=1 ) ! [IN]
       call FILEIO_read( TGL_URB(:,:,:),                                       & ! [OUT]
                         URBAN_RESTART_IN_BASENAME, 'TGL_URB', 'Urban', step=1 ) ! [IN]

       call URBAN_vars_fillhalo

       call URBAN_vars_total
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** restart file for urban is not specified.'
    endif

    return
  end subroutine URBAN_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Write urban restart
  subroutine URBAN_vars_restart_write
    use scale_time, only: &
       TIME_gettimelabel
    use scale_fileio, only: &
       FILEIO_write
    implicit none

    character(len=15)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    if ( URBAN_RESTART_OUT_BASENAME /= '' ) then

       call TIME_gettimelabel( timelabel )
       write(basename,'(A,A,A)') trim(URBAN_RESTART_OUT_BASENAME), '_', trim(timelabel)

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Output restart file (URBAN) ***'
       if( IO_L ) write(IO_FID_LOG,*) '*** filename: ', trim(basename)

       call URBAN_vars_total

       call FILEIO_write( TR_URB(:,:),  basename,                                           URBAN_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_TR_URB), VAR_DESC(I_TR_URB), VAR_UNIT(I_TR_URB), 'XY', URBAN_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( TB_URB(:,:),  basename,                                           URBAN_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_TB_URB), VAR_DESC(I_TB_URB), VAR_UNIT(I_TB_URB), 'XY', URBAN_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( TG_URB(:,:),  basename,                                           URBAN_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_TG_URB), VAR_DESC(I_TG_URB), VAR_UNIT(I_TG_URB), 'XY', URBAN_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( TC_URB(:,:),  basename,                                           URBAN_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_TC_URB), VAR_DESC(I_TC_URB), VAR_UNIT(I_TC_URB), 'XY', URBAN_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( QC_URB(:,:),  basename,                                           URBAN_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_QC_URB), VAR_DESC(I_QC_URB), VAR_UNIT(I_QC_URB), 'XY', URBAN_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( UC_URB(:,:),  basename,                                           URBAN_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_UC_URB), VAR_DESC(I_UC_URB), VAR_UNIT(I_UC_URB), 'XY', URBAN_RESTART_OUT_DTYPE  ) ! [IN]

       call FILEIO_write( TRL_URB(:,:,:), basename,                                               URBAN_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_TRL_URB), VAR_DESC(I_TRL_URB), VAR_UNIT(I_TRL_URB), 'Urban', URBAN_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( TBL_URB(:,:,:), basename,                                               URBAN_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_TBL_URB), VAR_DESC(I_TBL_URB), VAR_UNIT(I_TBL_URB), 'Urban', URBAN_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( TGL_URB(:,:,:), basename,                                               URBAN_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_TGL_URB), VAR_DESC(I_TGL_URB), VAR_UNIT(I_TGL_URB), 'Urban', URBAN_RESTART_OUT_DTYPE  ) ! [IN]

    endif

    return
  end subroutine URBAN_vars_restart_write

  !-----------------------------------------------------------------------------
  !> History output set for urban variables
  subroutine URBAN_vars_history
    use scale_time, only: &
       TIME_DTSEC_URBAN
    use scale_history, only: &
       HIST_in
    implicit none
    !---------------------------------------------------------------------------

    if ( URBAN_VARS_CHECKRANGE ) then
       call VALCHECK( TR_URB(:,:), 0.0_RP, 1000.0_RP, VAR_NAME(I_TR_URB), __FILE__, __LINE__ )
       call VALCHECK( TB_URB(:,:), 0.0_RP, 1000.0_RP, VAR_NAME(I_TB_URB), __FILE__, __LINE__ )
       call VALCHECK( TG_URB(:,:), 0.0_RP, 1000.0_RP, VAR_NAME(I_TG_URB), __FILE__, __LINE__ )
       call VALCHECK( TC_URB(:,:), 0.0_RP, 1000.0_RP, VAR_NAME(I_TC_URB), __FILE__, __LINE__ )
       call VALCHECK( QC_URB(:,:), 0.0_RP, 1000.0_RP, VAR_NAME(I_QC_URB), __FILE__, __LINE__ )
       call VALCHECK( UC_URB(:,:), 0.0_RP, 1000.0_RP, VAR_NAME(I_UC_URB), __FILE__, __LINE__ )

       call VALCHECK( TRL_URB(:,:,:), 0.0_RP, 1000.0_RP, VAR_NAME(I_TRL_URB), __FILE__, __LINE__ )
       call VALCHECK( TBL_URB(:,:,:), 0.0_RP, 1000.0_RP, VAR_NAME(I_TBL_URB), __FILE__, __LINE__ )
       call VALCHECK( TGL_URB(:,:,:), 0.0_RP, 1000.0_RP, VAR_NAME(I_TGL_URB), __FILE__, __LINE__ )
    endif

    call HIST_in( TR_URB(:,:), 'TR_URB', VAR_DESC(I_TR_URB), VAR_UNIT(I_TR_URB), TIME_DTSEC_URBAN )
    call HIST_in( TB_URB(:,:), 'TB_URB', VAR_DESC(I_TB_URB), VAR_UNIT(I_TB_URB), TIME_DTSEC_URBAN )
    call HIST_in( TG_URB(:,:), 'TG_URB', VAR_DESC(I_TG_URB), VAR_UNIT(I_TG_URB), TIME_DTSEC_URBAN )
    call HIST_in( TC_URB(:,:), 'TC_URB', VAR_DESC(I_TC_URB), VAR_UNIT(I_TC_URB), TIME_DTSEC_URBAN )
    call HIST_in( QC_URB(:,:), 'QC_URB', VAR_DESC(I_QC_URB), VAR_UNIT(I_QC_URB), TIME_DTSEC_URBAN )
    call HIST_in( UC_URB(:,:), 'UC_URB', VAR_DESC(I_UC_URB), VAR_UNIT(I_UC_URB), TIME_DTSEC_URBAN )

    call HIST_in( TRL_URB(:,:,:), 'TRL_URB', VAR_DESC(I_TRL_URB), VAR_UNIT(I_TRL_URB), TIME_DTSEC_URBAN, zdim='urban' )
    call HIST_in( TBL_URB(:,:,:), 'TBL_URB', VAR_DESC(I_TBL_URB), VAR_UNIT(I_TBL_URB), TIME_DTSEC_URBAN, zdim='urban' )
    call HIST_in( TGL_URB(:,:,:), 'TGL_URB', VAR_DESC(I_TGL_URB), VAR_UNIT(I_TGL_URB), TIME_DTSEC_URBAN, zdim='urban' )

    return
  end subroutine URBAN_vars_history

  !-----------------------------------------------------------------------------
  !> Budget monitor for urban
  subroutine URBAN_vars_total
    use scale_stats, only: &
       STAT_checktotal, &
       STAT_total
    implicit none

    real(RP) :: total
    integer  :: k
    !---------------------------------------------------------------------------

    if ( STAT_checktotal ) then
       call STAT_total( total, TR_URB(:,:), VAR_NAME(I_TR_URB) )
       call STAT_total( total, TB_URB(:,:), VAR_NAME(I_TB_URB) )
       call STAT_total( total, TG_URB(:,:), VAR_NAME(I_TG_URB) )
       call STAT_total( total, TC_URB(:,:), VAR_NAME(I_TC_URB) )
       call STAT_total( total, QC_URB(:,:), VAR_NAME(I_QC_URB) )
       call STAT_total( total, UC_URB(:,:), VAR_NAME(I_UC_URB) )

       do k = UKS, UKE
          call STAT_total( total, TRL_URB(k,:,:), VAR_NAME(I_TRL_URB) )
          call STAT_total( total, TBL_URB(k,:,:), VAR_NAME(I_TBL_URB) )
          call STAT_total( total, TGL_URB(k,:,:), VAR_NAME(I_TGL_URB) )
       enddo
    endif

    return
  end subroutine URBAN_vars_total

end module mod_urban_vars
