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
  public :: URBAN_vars_setup
  public :: URBAN_vars_fillhalo
  public :: URBAN_vars_restart_read
  public :: URBAN_vars_restart_write
  public :: URBAN_vars_history
  public :: URBAN_vars_total
  public :: URBAN_vars_external_in

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical, public :: URBAN_RESTART_OUTPUT = .false. !< output restart file?

  ! prognostic variables
  real(RP), public, allocatable :: TR_URB (:,:)   ! Surface temperature of roof [K]
  real(RP), public, allocatable :: TB_URB (:,:)   ! Surface temperature of wall [K]
  real(RP), public, allocatable :: TG_URB (:,:)   ! Surface temperature of road [K]
  real(RP), public, allocatable :: TC_URB (:,:)   ! Diagnostic canopy air temperature [K]
  real(RP), public, allocatable :: QC_URB (:,:)   ! Diagnostic canopy humidity [-]
  real(RP), public, allocatable :: UC_URB (:,:)   ! Diagnostic canopy wind [m/s]
  real(RP), public, allocatable :: TS_URB (:,:)   ! Diagnostic surface temperature [K]
  real(RP), public, allocatable :: SHR_URB (:,:)  ! Sensible heat flux from roof [W/m2]
  real(RP), public, allocatable :: SHB_URB (:,:)  ! Sensible heat flux from wall [W/m2]
  real(RP), public, allocatable :: SHG_URB (:,:)  ! Sensible heat flux from road [W/m2]
  real(RP), public, allocatable :: LHR_URB (:,:)  ! Latent heat flux from roof [W/m2]
  real(RP), public, allocatable :: LHB_URB (:,:)  ! Latent heat flux from wall [W/m2]
  real(RP), public, allocatable :: LHG_URB (:,:)  ! Latent heat flux from road [W/m2]
  real(RP), public, allocatable :: GHR_URB (:,:)  ! Ground heat flux on roof [W/m2]
  real(RP), public, allocatable :: GHB_URB (:,:)  ! Ground heat flux on wall [W/m2]
  real(RP), public, allocatable :: GHG_URB (:,:)  ! Ground heat flux on road [W/m2]
  real(RP), public, allocatable :: RnR_URB (:,:)  ! Net radiation on roof [W/m2]
  real(RP), public, allocatable :: RnB_URB (:,:)  ! Net radiation on wall [W/m2]
  real(RP), public, allocatable :: RnG_URB (:,:)  ! Net radiation on road [W/m2]
  real(RP), public, allocatable :: TRL_URB(:,:,:) ! temperature in layer of roof [K]
  real(RP), public, allocatable :: TBL_URB(:,:,:) ! temperature in layer of wall [K]
  real(RP), public, allocatable :: TGL_URB(:,:,:) ! temperature in layer of road [K]
  real(RP), public, allocatable :: RAINR_URB(:,:) ! Rain storage on roof [mm=kg/m2]
  real(RP), public, allocatable :: RAINB_URB(:,:) ! Rain storage on roof [mm=kg/m2]
  real(RP), public, allocatable :: RAING_URB(:,:) ! Rain storage on roof [mm=kg/m2]
  real(RP), public, allocatable :: ROFF_URB(:,:)  ! Runoff from urban [mm=kg/m2]


  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_LONG),  private :: URBAN_RESTART_IN_BASENAME  = ''              !< basename of the restart file
  character(len=H_LONG),  private :: URBAN_RESTART_OUT_BASENAME = ''              !< basename of the output file
  character(len=H_MID),   private :: URBAN_RESTART_OUT_TITLE    = 'URBAN restart' !< title    of the output file
  character(len=H_MID),   private :: URBAN_RESTART_OUT_DTYPE    = 'DEFAULT'       !< REAL4 or REAL8

  logical,                private :: URBAN_VARS_CHECKRANGE      = .false.

  integer,                private, parameter :: VMAX = 26
  integer,                private, parameter :: I_TR_URB  = 1
  integer,                private, parameter :: I_TB_URB  = 2
  integer,                private, parameter :: I_TG_URB  = 3
  integer,                private, parameter :: I_TC_URB  = 4
  integer,                private, parameter :: I_QC_URB  = 5
  integer,                private, parameter :: I_UC_URB  = 6
  integer,                private, parameter :: I_TS_URB  = 7
  integer,                private, parameter :: I_SHR_URB = 8
  integer,                private, parameter :: I_SHB_URB = 9
  integer,                private, parameter :: I_SHG_URB = 10
  integer,                private, parameter :: I_LHR_URB = 11
  integer,                private, parameter :: I_LHB_URB = 12
  integer,                private, parameter :: I_LHG_URB = 13
  integer,                private, parameter :: I_GHR_URB = 14
  integer,                private, parameter :: I_GHB_URB = 15
  integer,                private, parameter :: I_GHG_URB = 16
  integer,                private, parameter :: I_RnR_URB = 17
  integer,                private, parameter :: I_RnB_URB = 18
  integer,                private, parameter :: I_RnG_URB = 19
  integer,                private, parameter :: I_TRL_URB = 20
  integer,                private, parameter :: I_TBL_URB = 21
  integer,                private, parameter :: I_TGL_URB = 22
  integer,                private, parameter :: I_RAINR_URB = 23
  integer,                private, parameter :: I_RAINB_URB = 24
  integer,                private, parameter :: I_RAING_URB = 25
  integer,                private, parameter :: I_ROFF_URB  = 26

  character(len=H_SHORT), private            :: VAR_NAME(VMAX) !< name  of the urban variables
  character(len=H_MID),   private            :: VAR_DESC(VMAX) !< desc. of the urban variables
  character(len=H_SHORT), private            :: VAR_UNIT(VMAX) !< unit  of the urban variables

  data VAR_NAME / 'TR_URB' ,  &
                  'TB_URB' ,  &
                  'TG_URB' ,  &
                  'TC_URB' ,  &
                  'QC_URB' ,  &
                  'UC_URB' ,  &
                  'TS_URB' ,  &
                  'SHR_URB' , &
                  'SHB_URB' , &
                  'SHG_URB' , &
                  'LHR_URB' , &
                  'LHB_URB' , &
                  'LHG_URB' , &
                  'GHR_URB' , &
                  'GHB_URB' , &
                  'GHG_URB' , &
                  'RnR_URB' , &
                  'RnB_URB' , &
                  'RnG_URB' , &
                  'TRL_URB' , &
                  'TBL_URB' , &
                  'TGL_URB' , &
                  'RAINR_URB' , &
                  'RAINB_URB' , &
                  'RAING_URB' , &
                  'ROFF_URB'    /

  data VAR_DESC / 'Surface temperature of roof',       &
                  'Surface temperature of wall',       &
                  'Surface temperature of road',       &
                  'Diagnostic canopy air temperature', &
                  'Diagnostic canopy humidity',        &
                  'Diagnostic canopy wind',            &
                  'Diagnostic surface temperature',    &
                  'Sensible heat flux from roof',      &
                  'Sensible heat flux from wall',      &
                  'Sensible heat flux from road',      &
                  'Latent heat flux from roof',        &
                  'Latent heat flux from wall',        &
                  'Latent heat flux from road',        &
                  'Ground heat flux on roof',          &
                  'Ground heat flux on wall',          &
                  'Ground heat flux on road',          &
                  'Net radiation on roof',             &
                  'Net radiation on wall',             &
                  'Net radiation on road',             &
                  'temperature in layer of roof',      &
                  'temperature in layer of wall',      &
                  'temperature in layer of road',      &
                  'rain strage on roof',               &
                  'rain strage on building',           &
                  'rain strage on road',               &
                  'runoff from urban'                  /

  data VAR_UNIT / 'K',     &
                  'K',     &
                  'K',     &
                  'K',     &
                  '-',     &
                  'm/s',   &
                  'K',     &
                  'W/m2',  &
                  'W/m2',  &
                  'W/m2',  &
                  'W/m2',  &
                  'W/m2',  &
                  'W/m2',  &
                  'W/m2',  &
                  'W/m2',  &
                  'W/m2',  &
                  'W/m2',  &
                  'W/m2',  &
                  'W/m2',  &
                  'K',     &
                  'K',     &
                  'K',     &
                  'kg/m2', &
                  'kg/m2', &
                  'kg/m2', &
                  'kg/m2'  /

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
    integer :: iv
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
    allocate( TS_URB(IA,JA) )
    TR_URB (:,:) = UNDEF
    TB_URB (:,:) = UNDEF
    TG_URB (:,:) = UNDEF
    TC_URB (:,:) = UNDEF
    QC_URB (:,:) = UNDEF
    UC_URB (:,:) = UNDEF
    TS_URB (:,:) = UNDEF

    allocate( SHR_URB(IA,JA) )
    allocate( SHB_URB(IA,JA) )
    allocate( SHG_URB(IA,JA) )
    allocate( LHR_URB(IA,JA) )
    allocate( LHB_URB(IA,JA) )
    allocate( LHG_URB(IA,JA) )
    allocate( GHR_URB(IA,JA) )
    allocate( GHB_URB(IA,JA) )
    allocate( GHG_URB(IA,JA) )
    allocate( RnR_URB(IA,JA) )
    allocate( RnB_URB(IA,JA) )
    allocate( RnG_URB(IA,JA) )
    SHR_URB (:,:) = UNDEF
    SHB_URB (:,:) = UNDEF
    SHG_URB (:,:) = UNDEF
    LHR_URB (:,:) = UNDEF
    LHB_URB (:,:) = UNDEF
    LHG_URB (:,:) = UNDEF
    GHR_URB (:,:) = UNDEF
    GHB_URB (:,:) = UNDEF
    GHG_URB (:,:) = UNDEF
    RnR_URB (:,:) = UNDEF
    RnB_URB (:,:) = UNDEF
    RnG_URB (:,:) = UNDEF

    allocate( TRL_URB(UKS:UKE,IA,JA) )
    allocate( TBL_URB(UKS:UKE,IA,JA) )
    allocate( TGL_URB(UKS:UKE,IA,JA) )
    TRL_URB(:,:,:) = UNDEF
    TBL_URB(:,:,:) = UNDEF
    TGL_URB(:,:,:) = UNDEF

    allocate( RAINR_URB(IA,JA) )
    allocate( RAINB_URB(IA,JA) )
    allocate( RAING_URB(IA,JA) )
    allocate( ROFF_URB(IA,JA) )
    RAINR_URB (:,:) = 0.0_RP
    RAINB_URB (:,:) = 0.0_RP
    RAING_URB (:,:) = 0.0_RP
    ROFF_URB (:,:)  = 0.0_RP

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
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A15,A,A32,3(A))') &
               '***       |','VARNAME        ','|', 'DESCRIPTION                     ','[', 'UNIT            ',']'
    do iv = 1, VMAX
       if( IO_L ) write(IO_FID_LOG,'(1x,A,i3,A,A15,A,A32,3(A))') &
                  '*** NO.',iv,'|',VAR_NAME(iv),'|',VAR_DESC(iv),'[',VAR_UNIT(iv),']'
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
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart output? : NO'
       URBAN_RESTART_OUTPUT = .false.
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

    call COMM_vars8( TR_URB(:,:),    1+UKE*3 )
    call COMM_vars8( TB_URB(:,:),    2+UKE*3 )
    call COMM_vars8( TG_URB(:,:),    3+UKE*3 )
    call COMM_vars8( TC_URB(:,:),    4+UKE*3 )
    call COMM_vars8( QC_URB(:,:),    5+UKE*3 )
    call COMM_vars8( UC_URB(:,:),    6+UKE*3 )
    call COMM_vars8( TS_URB(:,:),    7+UKE*3 )
    call COMM_vars8( SHR_URB(:,:),   8+UKE*3 )
    call COMM_vars8( SHB_URB(:,:),   9+UKE*3 )
    call COMM_vars8( SHG_URB(:,:),   10+UKE*3 )
    call COMM_vars8( LHR_URB(:,:),   11+UKE*3 )
    call COMM_vars8( LHB_URB(:,:),   12+UKE*3 )
    call COMM_vars8( LHG_URB(:,:),   13+UKE*3 )
    call COMM_vars8( GHR_URB(:,:),   14+UKE*3 )
    call COMM_vars8( GHB_URB(:,:),   15+UKE*3 )
    call COMM_vars8( GHG_URB(:,:),   16+UKE*3 )
    call COMM_vars8( RnR_URB(:,:),   17+UKE*3 )
    call COMM_vars8( RnB_URB(:,:),   18+UKE*3 )
    call COMM_vars8( RnG_URB(:,:),   19+UKE*3 )
    call COMM_vars8( RAINR_URB(:,:), 20+UKE*3 )
    call COMM_vars8( RAINB_URB(:,:), 21+UKE*3 )
    call COMM_vars8( RAING_URB(:,:), 22+UKE*3 )
    call COMM_vars8( ROFF_URB(:,:),  23+UKE*3 )

    do k = UKS, UKE
      call COMM_wait ( tmp1(:,:,k), k       )
      call COMM_wait ( tmp2(:,:,k), k+UKE   )
      call COMM_wait ( tmp3(:,:,k), k+UKE*2 )
    end do

    call COMM_wait ( TR_URB(:,:),    1+UKE*3 )
    call COMM_wait ( TB_URB(:,:),    2+UKE*3 )
    call COMM_wait ( TG_URB(:,:),    3+UKE*3 )
    call COMM_wait ( TC_URB(:,:),    4+UKE*3 )
    call COMM_wait ( QC_URB(:,:),    5+UKE*3 )
    call COMM_wait ( UC_URB(:,:),    6+UKE*3 )
    call COMM_wait ( TS_URB(:,:),    7+UKE*3 )
    call COMM_wait ( SHR_URB(:,:),   8+UKE*3 )
    call COMM_wait ( SHB_URB(:,:),   9+UKE*3 )
    call COMM_wait ( SHG_URB(:,:),   10+UKE*3 )
    call COMM_wait ( LHR_URB(:,:),   11+UKE*3 )
    call COMM_wait ( LHB_URB(:,:),   12+UKE*3 )
    call COMM_wait ( LHG_URB(:,:),   13+UKE*3 )
    call COMM_wait ( GHR_URB(:,:),   14+UKE*3 )
    call COMM_wait ( GHB_URB(:,:),   15+UKE*3 )
    call COMM_wait ( GHG_URB(:,:),   16+UKE*3 )
    call COMM_wait ( RnR_URB(:,:),   17+UKE*3 )
    call COMM_wait ( RnB_URB(:,:),   18+UKE*3 )
    call COMM_wait ( RnG_URB(:,:),   19+UKE*3 )
    call COMM_wait ( RAINR_URB(:,:), 20+UKE*3 )
    call COMM_wait ( RAINB_URB(:,:), 21+UKE*3 )
    call COMM_wait ( RAING_URB(:,:), 22+UKE*3 )
    call COMM_wait ( ROFF_URB(:,:),  23+UKE*3 )

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
       call FILEIO_read( TS_URB(:,:),                                      & ! [OUT]
                         URBAN_RESTART_IN_BASENAME, 'TS_URB', 'XY', step=1 ) ! [IN]

       call FILEIO_read( TRL_URB(:,:,:),                                       & ! [OUT]
                         URBAN_RESTART_IN_BASENAME, 'TRL_URB', 'Urban', step=1 ) ! [IN]
       call FILEIO_read( TBL_URB(:,:,:),                                       & ! [OUT]
                         URBAN_RESTART_IN_BASENAME, 'TBL_URB', 'Urban', step=1 ) ! [IN]
       call FILEIO_read( TGL_URB(:,:,:),                                       & ! [OUT]
                         URBAN_RESTART_IN_BASENAME, 'TGL_URB', 'Urban', step=1 ) ! [IN]

       call FILEIO_read( RAINR_URB(:,:),                                      & ! [OUT]
                         URBAN_RESTART_IN_BASENAME, 'RAINR_URB', 'XY', step=1 ) ! [IN]
       call FILEIO_read( RAINB_URB(:,:),                                      & ! [OUT]
                         URBAN_RESTART_IN_BASENAME, 'RAINB_URB', 'XY', step=1 ) ! [IN]
       call FILEIO_read( RAING_URB(:,:),                                      & ! [OUT]
                         URBAN_RESTART_IN_BASENAME, 'RAING_URB', 'XY', step=1 ) ! [IN]
       call FILEIO_read( ROFF_URB(:,:),                                       & ! [OUT]
                         URBAN_RESTART_IN_BASENAME, 'ROFF_URB', 'XY', step=1 )  ! [IN]

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
       call FILEIO_write( TS_URB(:,:),  basename,                                           URBAN_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_TS_URB), VAR_DESC(I_TS_URB), VAR_UNIT(I_TS_URB), 'XY', URBAN_RESTART_OUT_DTYPE  ) ! [IN]

       call FILEIO_write( TRL_URB(:,:,:), basename,                                               URBAN_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_TRL_URB), VAR_DESC(I_TRL_URB), VAR_UNIT(I_TRL_URB), 'Urban', URBAN_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( TBL_URB(:,:,:), basename,                                               URBAN_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_TBL_URB), VAR_DESC(I_TBL_URB), VAR_UNIT(I_TBL_URB), 'Urban', URBAN_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( TGL_URB(:,:,:), basename,                                               URBAN_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_TGL_URB), VAR_DESC(I_TGL_URB), VAR_UNIT(I_TGL_URB), 'Urban', URBAN_RESTART_OUT_DTYPE  ) ! [IN]

       call FILEIO_write( RAINR_URB(:,:),  basename,                                                 URBAN_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_RAINR_URB), VAR_DESC(I_RAINR_URB), VAR_UNIT(I_RAINR_URB), 'XY', URBAN_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( RAINB_URB(:,:),  basename,                                                 URBAN_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_RAINB_URB), VAR_DESC(I_RAINB_URB), VAR_UNIT(I_RAINB_URB), 'XY', URBAN_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( RAING_URB(:,:),  basename,                                                 URBAN_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_RAING_URB), VAR_DESC(I_RAING_URB), VAR_UNIT(I_RAING_URB), 'XY', URBAN_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( ROFF_URB(:,:),  basename,                                               URBAN_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_ROFF_URB), VAR_DESC(I_ROFF_URB), VAR_UNIT(I_ROFF_URB), 'XY', URBAN_RESTART_OUT_DTYPE  ) ! [IN]

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
       call VALCHECK( RAINR_URB(:,:), -500.0_RP, 1000.0_RP, VAR_NAME(I_RAINR_URB), __FILE__, __LINE__ )
       call VALCHECK( RAINB_URB(:,:), -500.0_RP, 1000.0_RP, VAR_NAME(I_RAINB_URB), __FILE__, __LINE__ )
       call VALCHECK( RAING_URB(:,:), -500.0_RP, 1000.0_RP, VAR_NAME(I_RAING_URB), __FILE__, __LINE__ )
       call VALCHECK( ROFF_URB(:,:),  -500.0_RP, 1000.0_RP, VAR_NAME(I_ROFF_URB), __FILE__, __LINE__ )
    endif

    call HIST_in( TR_URB(:,:), 'TR_URB', VAR_DESC(I_TR_URB), VAR_UNIT(I_TR_URB), TIME_DTSEC_URBAN )
    call HIST_in( TB_URB(:,:), 'TB_URB', VAR_DESC(I_TB_URB), VAR_UNIT(I_TB_URB), TIME_DTSEC_URBAN )
    call HIST_in( TG_URB(:,:), 'TG_URB', VAR_DESC(I_TG_URB), VAR_UNIT(I_TG_URB), TIME_DTSEC_URBAN )
    call HIST_in( TC_URB(:,:), 'TC_URB', VAR_DESC(I_TC_URB), VAR_UNIT(I_TC_URB), TIME_DTSEC_URBAN )
    call HIST_in( QC_URB(:,:), 'QC_URB', VAR_DESC(I_QC_URB), VAR_UNIT(I_QC_URB), TIME_DTSEC_URBAN )
    call HIST_in( UC_URB(:,:), 'UC_URB', VAR_DESC(I_UC_URB), VAR_UNIT(I_UC_URB), TIME_DTSEC_URBAN )
    call HIST_in( TS_URB(:,:), 'TS_URB', VAR_DESC(I_TS_URB), VAR_UNIT(I_TS_URB), TIME_DTSEC_URBAN )

    call HIST_in( SHR_URB(:,:), 'SHR_URB', VAR_DESC(I_SHR_URB), VAR_UNIT(I_SHR_URB), TIME_DTSEC_URBAN )
    call HIST_in( SHB_URB(:,:), 'SHB_URB', VAR_DESC(I_SHB_URB), VAR_UNIT(I_SHB_URB), TIME_DTSEC_URBAN )
    call HIST_in( SHG_URB(:,:), 'SHG_URB', VAR_DESC(I_SHG_URB), VAR_UNIT(I_SHG_URB), TIME_DTSEC_URBAN )
    call HIST_in( LHR_URB(:,:), 'LHR_URB', VAR_DESC(I_LHR_URB), VAR_UNIT(I_LHR_URB), TIME_DTSEC_URBAN )
    call HIST_in( LHB_URB(:,:), 'LHB_URB', VAR_DESC(I_LHB_URB), VAR_UNIT(I_LHB_URB), TIME_DTSEC_URBAN )
    call HIST_in( LHG_URB(:,:), 'LHG_URB', VAR_DESC(I_LHG_URB), VAR_UNIT(I_LHG_URB), TIME_DTSEC_URBAN )
    call HIST_in( GHR_URB(:,:), 'GHR_URB', VAR_DESC(I_GHR_URB), VAR_UNIT(I_GHR_URB), TIME_DTSEC_URBAN )
    call HIST_in( GHB_URB(:,:), 'GHB_URB', VAR_DESC(I_GHB_URB), VAR_UNIT(I_GHB_URB), TIME_DTSEC_URBAN )
    call HIST_in( GHG_URB(:,:), 'GHG_URB', VAR_DESC(I_GHG_URB), VAR_UNIT(I_GHG_URB), TIME_DTSEC_URBAN )
    call HIST_in( RnR_URB(:,:), 'RnR_URB', VAR_DESC(I_RnR_URB), VAR_UNIT(I_RnR_URB), TIME_DTSEC_URBAN )
    call HIST_in( RnB_URB(:,:), 'RnB_URB', VAR_DESC(I_RnB_URB), VAR_UNIT(I_RnB_URB), TIME_DTSEC_URBAN )
    call HIST_in( RnG_URB(:,:), 'RnG_URB', VAR_DESC(I_RnG_URB), VAR_UNIT(I_RnG_URB), TIME_DTSEC_URBAN )

    call HIST_in( TRL_URB(:,:,:), 'TRL_URB', VAR_DESC(I_TRL_URB), VAR_UNIT(I_TRL_URB), TIME_DTSEC_URBAN, zdim='urban' )
    call HIST_in( TBL_URB(:,:,:), 'TBL_URB', VAR_DESC(I_TBL_URB), VAR_UNIT(I_TBL_URB), TIME_DTSEC_URBAN, zdim='urban' )
    call HIST_in( TGL_URB(:,:,:), 'TGL_URB', VAR_DESC(I_TGL_URB), VAR_UNIT(I_TGL_URB), TIME_DTSEC_URBAN, zdim='urban' )

    call HIST_in( RAINR_URB(:,:), 'RAINR_URB', VAR_DESC(I_RAINR_URB), VAR_UNIT(I_RAINR_URB), TIME_DTSEC_URBAN )
    call HIST_in( RAINB_URB(:,:), 'RAINB_URB', VAR_DESC(I_RAINB_URB), VAR_UNIT(I_RAINB_URB), TIME_DTSEC_URBAN )
    call HIST_in( RAING_URB(:,:), 'RAING_URB', VAR_DESC(I_RAING_URB), VAR_UNIT(I_RAING_URB), TIME_DTSEC_URBAN )
    call HIST_in( ROFF_URB(:,:),  'ROFF_URB',  VAR_DESC(I_ROFF_URB),  VAR_UNIT(I_ROFF_URB),  TIME_DTSEC_URBAN )

    return
  end subroutine URBAN_vars_history

  !-----------------------------------------------------------------------------
  !> Budget monitor for urban
  subroutine URBAN_vars_total
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STAT_total
    implicit none

    real(RP) :: total
    integer  :: k
    !---------------------------------------------------------------------------

    if ( STATISTICS_checktotal ) then
       call STAT_total( total, TR_URB(:,:), VAR_NAME(I_TR_URB) )
       call STAT_total( total, TB_URB(:,:), VAR_NAME(I_TB_URB) )
       call STAT_total( total, TG_URB(:,:), VAR_NAME(I_TG_URB) )
       call STAT_total( total, TC_URB(:,:), VAR_NAME(I_TC_URB) )
       call STAT_total( total, QC_URB(:,:), VAR_NAME(I_QC_URB) )
       call STAT_total( total, UC_URB(:,:), VAR_NAME(I_UC_URB) )
       call STAT_total( total, TS_URB(:,:), VAR_NAME(I_TS_URB) )
       call STAT_total( total, SHR_URB(:,:), VAR_NAME(I_SHR_URB) )
       call STAT_total( total, SHB_URB(:,:), VAR_NAME(I_SHB_URB) )
       call STAT_total( total, SHG_URB(:,:), VAR_NAME(I_SHG_URB) )
       call STAT_total( total, LHR_URB(:,:), VAR_NAME(I_LHR_URB) )
       call STAT_total( total, LHB_URB(:,:), VAR_NAME(I_LHB_URB) )
       call STAT_total( total, LHG_URB(:,:), VAR_NAME(I_LHG_URB) )
       call STAT_total( total, GHR_URB(:,:), VAR_NAME(I_GHR_URB) )
       call STAT_total( total, GHB_URB(:,:), VAR_NAME(I_GHB_URB) )
       call STAT_total( total, GHG_URB(:,:), VAR_NAME(I_GHG_URB) )
       call STAT_total( total, RnR_URB(:,:), VAR_NAME(I_RnR_URB) )
       call STAT_total( total, RnB_URB(:,:), VAR_NAME(I_RnB_URB) )
       call STAT_total( total, RnG_URB(:,:), VAR_NAME(I_RnG_URB) )

       do k = UKS, UKE
          call STAT_total( total, TRL_URB(k,:,:), VAR_NAME(I_TRL_URB) )
          call STAT_total( total, TBL_URB(k,:,:), VAR_NAME(I_TBL_URB) )
          call STAT_total( total, TGL_URB(k,:,:), VAR_NAME(I_TGL_URB) )
       enddo

       call STAT_total( total, RAINR_URB(:,:), VAR_NAME(I_RAINR_URB) )
       call STAT_total( total, RAINB_URB(:,:), VAR_NAME(I_RAINB_URB) )
       call STAT_total( total, RAING_URB(:,:), VAR_NAME(I_RAING_URB) )
       call STAT_total( total, ROFF_URB(:,:),  VAR_NAME(I_ROFF_URB) )
    endif

    return
  end subroutine URBAN_vars_total

  !-----------------------------------------------------------------------------
  !> Input from External I/O
  subroutine URBAN_vars_external_in( &
       ts_urb_in )
    implicit none

    real(RP), intent(in) :: ts_urb_in(IA,JA)
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** External Input (coupler) ***'

    TS_URB(:,:) = ts_urb_in(:,:)

    call URBAN_vars_fillhalo

    call URBAN_vars_total

    return
  end subroutine URBAN_vars_external_in

end module mod_urban_vars
