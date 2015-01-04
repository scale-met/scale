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
  public :: URBAN_vars_restart_read
  public :: URBAN_vars_restart_write
  public :: URBAN_vars_history
  public :: URBAN_vars_total
  public :: URBAN_vars_external_in

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical,               public :: URBAN_RESTART_OUTPUT       = .false.         !< output restart file?

  character(len=H_LONG), public :: URBAN_RESTART_IN_BASENAME  = ''              !< basename of the restart file
  character(len=H_LONG), public :: URBAN_RESTART_OUT_BASENAME = ''              !< basename of the output file
  character(len=H_MID),  public :: URBAN_RESTART_OUT_TITLE    = 'URBAN restart' !< title    of the output file
  character(len=H_MID),  public :: URBAN_RESTART_OUT_DTYPE    = 'DEFAULT'       !< REAL4 or REAL8

  ! prognostic variables
  real(RP), public, allocatable :: URBAN_TR   (:,:)   ! urban surface temperature of roof [K]
  real(RP), public, allocatable :: URBAN_TB   (:,:)   ! urban surface temperature of wall [K]
  real(RP), public, allocatable :: URBAN_TG   (:,:)   ! urban surface temperature of road [K]
  real(RP), public, allocatable :: URBAN_TC   (:,:)   ! urban canopy air temperature [K]
  real(RP), public, allocatable :: URBAN_QC   (:,:)   ! urban canopy humidity [kg/kg]
  real(RP), public, allocatable :: URBAN_UC   (:,:)   ! urban canopy wind [m/s]
  real(RP), public, allocatable :: URBAN_TRL  (:,:,:) ! urban temperature in layer of roof [K]
  real(RP), public, allocatable :: URBAN_TBL  (:,:,:) ! urban temperature in layer of wall [K]
  real(RP), public, allocatable :: URBAN_TGL  (:,:,:) ! urban temperature in layer of road [K]
  real(RP), public, allocatable :: URBAN_RAINR(:,:)   ! urban rain storage on roof [mm=kg/m2]
  real(RP), public, allocatable :: URBAN_RAINB(:,:)   ! urban rain storage on wall [mm=kg/m2]
  real(RP), public, allocatable :: URBAN_RAING(:,:)   ! urban rain storage on road [mm=kg/m2]
  real(RP), public, allocatable :: URBAN_ROFF (:,:)   ! urban runoff [mm=kg/m2]

  ! tendency variables
  real(RP), public, allocatable :: URBAN_TRL_t  (:,:,:) ! tendency of URBAN_TRL
  real(RP), public, allocatable :: URBAN_TBL_t  (:,:,:) ! tendency of URBAN_TBL
  real(RP), public, allocatable :: URBAN_TGL_t  (:,:,:) ! tendency of URBAN_TGL
  real(RP), public, allocatable :: URBAN_TC_t   (:,:)   ! tendency of URBAN_TC
  real(RP), public, allocatable :: URBAN_UC_t   (:,:)   ! tendency of URBAN_UC
  real(RP), public, allocatable :: URBAN_QC_t   (:,:)   ! tendency of URBAN_QC
  real(RP), public, allocatable :: URBAN_TR_t   (:,:)   ! tendency of URBAN_TR
  real(RP), public, allocatable :: URBAN_TB_t   (:,:)   ! tendency of URBAN_TB
  real(RP), public, allocatable :: URBAN_TG_t   (:,:)   ! tendency of URBAN_TG
  real(RP), public, allocatable :: URBAN_RAINR_t(:,:)   ! tendency of URBAN_RAINR
  real(RP), public, allocatable :: URBAN_RAINB_t(:,:)   ! tendency of URBAN_RAINB
  real(RP), public, allocatable :: URBAN_RAING_t(:,:)   ! tendency of URBAN_RAING
  real(RP), public, allocatable :: URBAN_ROFF_t (:,:)   ! tendency of URBAN_ROFF

  ! for restart
  real(RP), public, allocatable :: URBAN_SFC_TEMP  (:,:)   ! urban grid average of surface temperature [K]
  real(RP), public, allocatable :: URBAN_SFC_albedo(:,:,:) ! urban grid average of albedo [0-1]
  real(RP), public, allocatable :: URBAN_SFLX_MW   (:,:)   ! urban grid average of w-momentum flux [kg/m2/s]
  real(RP), public, allocatable :: URBAN_SFLX_MU   (:,:)   ! urban grid average of u-momentum flux [kg/m2/s]
  real(RP), public, allocatable :: URBAN_SFLX_MV   (:,:)   ! urban grid average of v-momentum flux [kg/m2/s]
  real(RP), public, allocatable :: URBAN_SFLX_SH   (:,:)   ! urban grid average of sensible heat flux [W/m2]
  real(RP), public, allocatable :: URBAN_SFLX_LH   (:,:)   ! urban grid average of latent heat flux [W/m2]
  real(RP), public, allocatable :: URBAN_SFLX_GH   (:,:)   ! urban grid average of ground heat flux [W/m2]
  real(RP), public, allocatable :: URBAN_SFLX_evap (:,:)   ! urban grid average of water vapor flux [kg/m2/s]

  ! diagnostic variables
  real(RP), public, allocatable :: URBAN_Z0M(:,:) ! urban grid average of rougness length (momentum) [m]
  real(RP), public, allocatable :: URBAN_Z0H(:,:) ! urban grid average of rougness length (heat) [m]
  real(RP), public, allocatable :: URBAN_Z0E(:,:) ! urban grid average of rougness length (vapor) [m]
  real(RP), public, allocatable :: URBAN_U10(:,:) ! urban grid average of velocity u at 10m [m/s]
  real(RP), public, allocatable :: URBAN_V10(:,:) ! urban grid average of velocity v at 10m [m/s]
  real(RP), public, allocatable :: URBAN_T2 (:,:) ! urban grid average of temperature at 2m [K]
  real(RP), public, allocatable :: URBAN_Q2 (:,:) ! urban grid average of water vapor at 2m [kg/kg]

  ! recieved atmospheric variables
  real(RP), public, allocatable :: ATMOS_TEMP     (:,:)
  real(RP), public, allocatable :: ATMOS_PRES     (:,:)
  real(RP), public, allocatable :: ATMOS_W        (:,:)
  real(RP), public, allocatable :: ATMOS_U        (:,:)
  real(RP), public, allocatable :: ATMOS_V        (:,:)
  real(RP), public, allocatable :: ATMOS_DENS     (:,:)
  real(RP), public, allocatable :: ATMOS_QV       (:,:)
  real(RP), public, allocatable :: ATMOS_PBL      (:,:)
  real(RP), public, allocatable :: ATMOS_SFC_PRES (:,:)
  real(RP), public, allocatable :: ATMOS_SFLX_LW  (:,:)
  real(RP), public, allocatable :: ATMOS_SFLX_SW  (:,:)
  real(RP), public, allocatable :: ATMOS_SFLX_prec(:,:)

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical,                private :: URBAN_VARS_CHECKRANGE      = .false.

  integer,                private, parameter :: VMAX        = 23
  integer,                private, parameter :: I_TR        = 1
  integer,                private, parameter :: I_TB        = 2
  integer,                private, parameter :: I_TG        = 3
  integer,                private, parameter :: I_TC        = 4
  integer,                private, parameter :: I_QC        = 5
  integer,                private, parameter :: I_UC        = 6
  integer,                private, parameter :: I_TRL       = 7
  integer,                private, parameter :: I_TBL       = 8
  integer,                private, parameter :: I_TGL       = 9
  integer,                private, parameter :: I_RAINR     = 10
  integer,                private, parameter :: I_RAINB     = 11
  integer,                private, parameter :: I_RAING     = 12
  integer,                private, parameter :: I_ROFF      = 13
  integer,                private, parameter :: I_SFC_TEMP  = 14
  integer,                private, parameter :: I_ALB_LW    = 15
  integer,                private, parameter :: I_ALB_SW    = 16
  integer,                private, parameter :: I_SFLX_MW   = 17
  integer,                private, parameter :: I_SFLX_MU   = 18
  integer,                private, parameter :: I_SFLX_MV   = 19
  integer,                private, parameter :: I_SFLX_SH   = 20
  integer,                private, parameter :: I_SFLX_LH   = 21
  integer,                private, parameter :: I_SFLX_GH   = 22
  integer,                private, parameter :: I_SFLX_evap = 23

  character(len=H_SHORT), private            :: VAR_NAME(VMAX) !< name  of the urban variables
  character(len=H_MID),   private            :: VAR_DESC(VMAX) !< desc. of the urban variables
  character(len=H_SHORT), private            :: VAR_UNIT(VMAX) !< unit  of the urban variables

  data VAR_NAME / 'URBAN_TR' ,       &
                  'URBAN_TB' ,       &
                  'URBAN_TG' ,       &
                  'URBAN_TC' ,       &
                  'URBAN_QC' ,       &
                  'URBAN_UC' ,       &
                  'URBAN_TRL' ,      &
                  'URBAN_TBL' ,      &
                  'URBAN_TGL' ,      &
                  'URBAN_RAINR' ,    &
                  'URBAN_RAINB' ,    &
                  'URBAN_RAING' ,    &
                  'URBAN_ROFF',      &
                  'URBAN_SFC_TEMP',  &
                  'URBAN_ALB_LW',    &
                  'URBAN_ALB_SW',    &
                  'URBAN_SFLX_MW',   &
                  'URBAN_SFLX_MU',   &
                  'URBAN_SFLX_MV',   &
                  'URBAN_SFLX_SH',   &
                  'URBAN_SFLX_LH',   &
                  'URBAN_SFLX_GH',   &
                  'URBAN_SFLX_evap'  /

  data VAR_DESC / 'urban surface temperature of roof',                &
                  'urban surface temperature of wall',                &
                  'urban surface temperature of road',                &
                  'urban canopy air temperature',                     &
                  'urban canopy humidity',                            &
                  'urban canopy wind',                                &
                  'urban temperature in layer of roof',               &
                  'urban temperature in layer of wall',               &
                  'urban temperature in layer of road',               &
                  'urban rain strage on roof',                        &
                  'urban rain strage on wall',                        &
                  'urban rain strage on road',                        &
                  'urban runoff ',                                    &
                  'urban grid average of temperature',                &
                  'urban grid average of albedo LW',                  &
                  'urban grid average of albedo SW',                  &
                  'urban grid average of w-momentum flux',            &
                  'urban grid average of u-momentum flux',            &
                  'urban grid average of v-momentum flux',            &
                  'urban grid average of sensible heat flux',         &
                  'urban grid average of latent heat flux',           &
                  'urban grid average of ground heat flux',           &
                  'urban grid average of water vapor flux'            /

  data VAR_UNIT / 'K',       &
                  'K',       &
                  'K',       &
                  'K',       &
                  'kg/kg',   &
                  'm/s',     &
                  'K',       &
                  'K',       &
                  'K',       &
                  'kg/m2',   &
                  'kg/m2',   &
                  'kg/m2',   &
                  'kg/m2',   &
                  'K',       &
                  '0-1',     &
                  '0-1',     &
                  'kg/m2/s', &
                  'kg/m2/s', &
                  'kg/m2/s', &
                  'W/m2',    &
                  'W/m2',    &
                  'W/m2',    &
                  'kg/m2/s'  /

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

    allocate( URBAN_TR   (IA,JA)         )
    allocate( URBAN_TB   (IA,JA)         )
    allocate( URBAN_TG   (IA,JA)         )
    allocate( URBAN_TC   (IA,JA)         )
    allocate( URBAN_QC   (IA,JA)         )
    allocate( URBAN_UC   (IA,JA)         )
    allocate( URBAN_TRL  (UKS:UKE,IA,JA) )
    allocate( URBAN_TBL  (UKS:UKE,IA,JA) )
    allocate( URBAN_TGL  (UKS:UKE,IA,JA) )
    allocate( URBAN_RAINR(IA,JA)         )
    allocate( URBAN_RAINB(IA,JA)         )
    allocate( URBAN_RAING(IA,JA)         )
    allocate( URBAN_ROFF (IA,JA)         )
    URBAN_TR   (:,:)   = UNDEF
    URBAN_TB   (:,:)   = UNDEF
    URBAN_TG   (:,:)   = UNDEF
    URBAN_TC   (:,:)   = UNDEF
    URBAN_QC   (:,:)   = UNDEF
    URBAN_UC   (:,:)   = UNDEF
    URBAN_TRL  (:,:,:) = UNDEF
    URBAN_TBL  (:,:,:) = UNDEF
    URBAN_TGL  (:,:,:) = UNDEF
    URBAN_RAINR(:,:)   = UNDEF
    URBAN_RAINB(:,:)   = UNDEF
    URBAN_RAING(:,:)   = UNDEF
    URBAN_ROFF (:,:)   = UNDEF

    allocate( URBAN_TR_t   (IA,JA)         )
    allocate( URBAN_TB_t   (IA,JA)         )
    allocate( URBAN_TG_t   (IA,JA)         )
    allocate( URBAN_TC_t   (IA,JA)         )
    allocate( URBAN_QC_t   (IA,JA)         )
    allocate( URBAN_UC_t   (IA,JA)         )
    allocate( URBAN_TRL_t  (UKS:UKE,IA,JA) )
    allocate( URBAN_TBL_t  (UKS:UKE,IA,JA) )
    allocate( URBAN_TGL_t  (UKS:UKE,IA,JA) )
    allocate( URBAN_RAINR_t(IA,JA)         )
    allocate( URBAN_RAINB_t(IA,JA)         )
    allocate( URBAN_RAING_t(IA,JA)         )
    allocate( URBAN_ROFF_t (IA,JA)         )
    URBAN_TR_t   (:,:)   = UNDEF
    URBAN_TB_t   (:,:)   = UNDEF
    URBAN_TG_t   (:,:)   = UNDEF
    URBAN_TC_t   (:,:)   = UNDEF
    URBAN_QC_t   (:,:)   = UNDEF
    URBAN_UC_t   (:,:)   = UNDEF
    URBAN_TRL_t  (:,:,:) = UNDEF
    URBAN_TBL_t  (:,:,:) = UNDEF
    URBAN_TGL_t  (:,:,:) = UNDEF
    URBAN_RAINR_t(:,:)   = UNDEF
    URBAN_RAINB_t(:,:)   = UNDEF
    URBAN_RAING_t(:,:)   = UNDEF
    URBAN_ROFF_t (:,:)   = UNDEF

    allocate( URBAN_SFC_TEMP  (IA,JA)   )
    allocate( URBAN_SFC_albedo(IA,JA,2) )
    allocate( URBAN_SFLX_MW   (IA,JA)   )
    allocate( URBAN_SFLX_MU   (IA,JA)   )
    allocate( URBAN_SFLX_MV   (IA,JA)   )
    allocate( URBAN_SFLX_SH   (IA,JA)   )
    allocate( URBAN_SFLX_LH   (IA,JA)   )
    allocate( URBAN_SFLX_GH   (IA,JA)   )
    allocate( URBAN_SFLX_evap (IA,JA)   )
    URBAN_SFC_TEMP  (:,:)   = UNDEF
    URBAN_SFC_albedo(:,:,:) = UNDEF
    URBAN_SFLX_MW   (:,:)   = UNDEF
    URBAN_SFLX_MU   (:,:)   = UNDEF
    URBAN_SFLX_MV   (:,:)   = UNDEF
    URBAN_SFLX_SH   (:,:)   = UNDEF
    URBAN_SFLX_LH   (:,:)   = UNDEF
    URBAN_SFLX_GH   (:,:)   = UNDEF
    URBAN_SFLX_evap (:,:)   = UNDEF

    allocate( URBAN_Z0M(IA,JA) )
    allocate( URBAN_Z0H(IA,JA) )
    allocate( URBAN_Z0E(IA,JA) )
    allocate( URBAN_U10(IA,JA) )
    allocate( URBAN_V10(IA,JA) )
    allocate( URBAN_T2 (IA,JA) )
    allocate( URBAN_Q2 (IA,JA) )
    URBAN_Z0M(:,:) = UNDEF
    URBAN_Z0H(:,:) = UNDEF
    URBAN_Z0E(:,:) = UNDEF
    URBAN_U10(:,:) = UNDEF
    URBAN_V10(:,:) = UNDEF
    URBAN_T2 (:,:) = UNDEF
    URBAN_Q2 (:,:) = UNDEF

    allocate( ATMOS_TEMP     (IA,JA) )
    allocate( ATMOS_PRES     (IA,JA) )
    allocate( ATMOS_W        (IA,JA) )
    allocate( ATMOS_U        (IA,JA) )
    allocate( ATMOS_V        (IA,JA) )
    allocate( ATMOS_DENS     (IA,JA) )
    allocate( ATMOS_QV       (IA,JA) )
    allocate( ATMOS_PBL      (IA,JA) )
    allocate( ATMOS_SFC_PRES (IA,JA) )
    allocate( ATMOS_SFLX_LW  (IA,JA) )
    allocate( ATMOS_SFLX_SW  (IA,JA) )
    allocate( ATMOS_SFLX_prec(IA,JA) )
    ATMOS_TEMP     (:,:) = UNDEF
    ATMOS_PRES     (:,:) = UNDEF
    ATMOS_W        (:,:) = UNDEF
    ATMOS_U        (:,:) = UNDEF
    ATMOS_V        (:,:) = UNDEF
    ATMOS_DENS     (:,:) = UNDEF
    ATMOS_QV       (:,:) = UNDEF
    ATMOS_PBL      (:,:) = UNDEF
    ATMOS_SFC_PRES (:,:) = UNDEF
    ATMOS_SFLX_LW  (:,:) = UNDEF
    ATMOS_SFLX_SW  (:,:) = UNDEF
    ATMOS_SFLX_prec(:,:) = UNDEF

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
  !> Read urban restart
  subroutine URBAN_vars_restart_read
    use scale_fileio, only: &
       FILEIO_read
    use mod_urban_admin, only: &
       URBAN_sw
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input restart file (URBAN) ***'

    if ( URBAN_sw .and. URBAN_RESTART_IN_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(URBAN_RESTART_IN_BASENAME)

       call FILEIO_read( URBAN_TR(:,:),                                          & ! [OUT]
                         URBAN_RESTART_IN_BASENAME, VAR_NAME(I_TR), 'XY', step=1 ) ! [IN]
       call FILEIO_read( URBAN_TB(:,:),                                          & ! [OUT]
                         URBAN_RESTART_IN_BASENAME, VAR_NAME(I_TB), 'XY', step=1 ) ! [IN]
       call FILEIO_read( URBAN_TG(:,:),                                          & ! [OUT]
                         URBAN_RESTART_IN_BASENAME, VAR_NAME(I_TG), 'XY', step=1 ) ! [IN]
       call FILEIO_read( URBAN_TC(:,:),                                          & ! [OUT]
                         URBAN_RESTART_IN_BASENAME, VAR_NAME(I_TC), 'XY', step=1 ) ! [IN]
       call FILEIO_read( URBAN_QC(:,:),                                          & ! [OUT]
                         URBAN_RESTART_IN_BASENAME, VAR_NAME(I_QC), 'XY', step=1 ) ! [IN]
       call FILEIO_read( URBAN_UC(:,:),                                          & ! [OUT]
                         URBAN_RESTART_IN_BASENAME, VAR_NAME(I_UC), 'XY', step=1 ) ! [IN]

       call FILEIO_read( URBAN_TRL(:,:,:),                                           & ! [OUT]
                         URBAN_RESTART_IN_BASENAME, VAR_NAME(I_TRL), 'Urban', step=1 ) ! [IN]
       call FILEIO_read( URBAN_TBL(:,:,:),                                           & ! [OUT]
                         URBAN_RESTART_IN_BASENAME, VAR_NAME(I_TBL), 'Urban', step=1 ) ! [IN]
       call FILEIO_read( URBAN_TGL(:,:,:),                                           & ! [OUT]
                         URBAN_RESTART_IN_BASENAME, VAR_NAME(I_TGL), 'Urban', step=1 ) ! [IN]

       call FILEIO_read( URBAN_RAINR(:,:),                                          & ! [OUT]
                         URBAN_RESTART_IN_BASENAME, VAR_NAME(I_RAINR), 'XY', step=1 ) ! [IN]
       call FILEIO_read( URBAN_RAINB(:,:),                                          & ! [OUT]
                         URBAN_RESTART_IN_BASENAME, VAR_NAME(I_RAINB), 'XY', step=1 ) ! [IN]
       call FILEIO_read( URBAN_RAING(:,:),                                          & ! [OUT]
                         URBAN_RESTART_IN_BASENAME, VAR_NAME(I_RAING), 'XY', step=1 ) ! [IN]
       call FILEIO_read( URBAN_ROFF(:,:),                                           & ! [OUT]
                         URBAN_RESTART_IN_BASENAME, VAR_NAME(I_ROFF),  'XY', step=1 ) ! [IN]

       call FILEIO_read( URBAN_SFC_TEMP(:,:),                                          & ! [OUT]
                         URBAN_RESTART_IN_BASENAME, VAR_NAME(I_SFC_TEMP), 'XY', step=1 ) ! [IN]
       call FILEIO_read( URBAN_SFC_albedo(:,:,I_LW),                                   & ! [OUT]
                         URBAN_RESTART_IN_BASENAME, VAR_NAME(I_ALB_LW),   'XY', step=1 ) ! [IN]
       call FILEIO_read( URBAN_SFC_albedo(:,:,I_SW),                                   & ! [OUT]
                         URBAN_RESTART_IN_BASENAME, VAR_NAME(I_ALB_SW),   'XY', step=1 ) ! [IN]

       call FILEIO_read( URBAN_SFLX_MW(:,:),                                            & ! [OUT]
                         URBAN_RESTART_IN_BASENAME, VAR_NAME(I_SFLX_MW),   'XY', step=1 ) ! [IN]
       call FILEIO_read( URBAN_SFLX_MU(:,:),                                            & ! [OUT]
                         URBAN_RESTART_IN_BASENAME, VAR_NAME(I_SFLX_MU),   'XY', step=1 ) ! [IN]
       call FILEIO_read( URBAN_SFLX_MV(:,:),                                            & ! [OUT]
                         URBAN_RESTART_IN_BASENAME, VAR_NAME(I_SFLX_MV),   'XY', step=1 ) ! [IN]
       call FILEIO_read( URBAN_SFLX_SH(:,:),                                            & ! [OUT]
                         URBAN_RESTART_IN_BASENAME, VAR_NAME(I_SFLX_SH),   'XY', step=1 ) ! [IN]
       call FILEIO_read( URBAN_SFLX_LH(:,:),                                            & ! [OUT]
                         URBAN_RESTART_IN_BASENAME, VAR_NAME(I_SFLX_LH),   'XY', step=1 ) ! [IN]
       call FILEIO_read( URBAN_SFLX_GH(:,:),                                            & ! [OUT]
                         URBAN_RESTART_IN_BASENAME, VAR_NAME(I_SFLX_GH),   'XY', step=1 ) ! [IN]
       call FILEIO_read( URBAN_SFLX_evap(:,:),                                          & ! [OUT]
                         URBAN_RESTART_IN_BASENAME, VAR_NAME(I_SFLX_evap), 'XY', step=1 ) ! [IN]

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
    use mod_urban_admin, only: &
       URBAN_sw
    implicit none

    character(len=15)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    if ( URBAN_sw .and. URBAN_RESTART_OUT_BASENAME /= '' ) then

       call TIME_gettimelabel( timelabel )
       write(basename,'(A,A,A)') trim(URBAN_RESTART_OUT_BASENAME), '_', trim(timelabel)

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Output restart file (URBAN) ***'
       if( IO_L ) write(IO_FID_LOG,*) '*** filename: ', trim(basename)

       call URBAN_vars_total

       call FILEIO_write( URBAN_TR(:,:),  basename, URBAN_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_TR), VAR_DESC(I_TR), VAR_UNIT(I_TR),    & ! [IN]
                          'XY', URBAN_RESTART_OUT_DTYPE, nohalo=.true.       ) ! [IN]
       call FILEIO_write( URBAN_TB(:,:),  basename, URBAN_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_TB), VAR_DESC(I_TB), VAR_UNIT(I_TB),    & ! [IN]
                          'XY', URBAN_RESTART_OUT_DTYPE, nohalo=.true.       ) ! [IN]
       call FILEIO_write( URBAN_TG(:,:),  basename, URBAN_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_TG), VAR_DESC(I_TG), VAR_UNIT(I_TG),    & ! [IN]
                          'XY', URBAN_RESTART_OUT_DTYPE, nohalo=.true.       ) ! [IN]
       call FILEIO_write( URBAN_TC(:,:),  basename, URBAN_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_TC), VAR_DESC(I_TC), VAR_UNIT(I_TC),    & ! [IN]
                          'XY', URBAN_RESTART_OUT_DTYPE, nohalo=.true.       ) ! [IN]
       call FILEIO_write( URBAN_QC(:,:),  basename, URBAN_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_QC), VAR_DESC(I_QC), VAR_UNIT(I_QC),    & ! [IN]
                          'XY', URBAN_RESTART_OUT_DTYPE, nohalo=.true.       ) ! [IN]
       call FILEIO_write( URBAN_UC(:,:),  basename, URBAN_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_UC), VAR_DESC(I_UC), VAR_UNIT(I_UC),    & ! [IN]
                          'XY', URBAN_RESTART_OUT_DTYPE, nohalo=.true.       ) ! [IN]

       call FILEIO_write( URBAN_TRL(:,:,:), basename, URBAN_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_TRL), VAR_DESC(I_TRL), VAR_UNIT(I_TRL),   & ! [IN]
                          'Urban', URBAN_RESTART_OUT_DTYPE, nohalo=.true.      ) ! [IN]
       call FILEIO_write( URBAN_TBL(:,:,:), basename, URBAN_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_TBL), VAR_DESC(I_TBL), VAR_UNIT(I_TBL),   & ! [IN]
                          'Urban', URBAN_RESTART_OUT_DTYPE, nohalo=.true.      ) ! [IN]
       call FILEIO_write( URBAN_TGL(:,:,:), basename, URBAN_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(I_TGL), VAR_DESC(I_TGL), VAR_UNIT(I_TGL),   & ! [IN]
                          'Urban', URBAN_RESTART_OUT_DTYPE, nohalo=.true.      ) ! [IN]

       call FILEIO_write( URBAN_RAINR(:,:),  basename, URBAN_RESTART_OUT_TITLE,    & ! [IN]
                          VAR_NAME(I_RAINR), VAR_DESC(I_RAINR), VAR_UNIT(I_RAINR), & ! [IN]
                          'XY', URBAN_RESTART_OUT_DTYPE, nohalo=.true.             ) ! [IN]
       call FILEIO_write( URBAN_RAINB(:,:),  basename, URBAN_RESTART_OUT_TITLE,    & ! [IN]
                          VAR_NAME(I_RAINB), VAR_DESC(I_RAINB), VAR_UNIT(I_RAINB), & ! [IN]
                          'XY', URBAN_RESTART_OUT_DTYPE, nohalo=.true.             ) ! [IN]
       call FILEIO_write( URBAN_RAING(:,:),  basename, URBAN_RESTART_OUT_TITLE,    & ! [IN]
                          VAR_NAME(I_RAING), VAR_DESC(I_RAING), VAR_UNIT(I_RAING), & ! [IN]
                          'XY', URBAN_RESTART_OUT_DTYPE, nohalo=.true.             ) ! [IN]
       call FILEIO_write( URBAN_ROFF(:,:),   basename, URBAN_RESTART_OUT_TITLE,    & ! [IN]
                          VAR_NAME(I_ROFF),  VAR_DESC(I_ROFF),  VAR_UNIT(I_ROFF),  & ! [IN]
                          'XY', URBAN_RESTART_OUT_DTYPE, nohalo=.true.             ) ! [IN]

       call FILEIO_write( URBAN_SFC_TEMP(:,:), basename, URBAN_RESTART_OUT_TITLE,           & ! [IN]
                          VAR_NAME(I_SFC_TEMP), VAR_DESC(I_SFC_TEMP), VAR_UNIT(I_SFC_TEMP), & ! [IN]
                          'XY', URBAN_RESTART_OUT_DTYPE, nohalo=.true.                      ) ! [IN]
       call FILEIO_write( URBAN_SFC_albedo(:,:,I_LW), basename, URBAN_RESTART_OUT_TITLE,    & ! [IN]
                          VAR_NAME(I_ALB_LW), VAR_DESC(I_ALB_LW), VAR_UNIT(I_ALB_LW),       & ! [IN]
                          'XY', URBAN_RESTART_OUT_DTYPE, nohalo=.true.                      ) ! [IN]
       call FILEIO_write( URBAN_SFC_albedo(:,:,I_SW), basename, URBAN_RESTART_OUT_TITLE,    & ! [IN]
                          VAR_NAME(I_ALB_SW), VAR_DESC(I_ALB_SW), VAR_UNIT(I_ALB_SW),       & ! [IN]
                          'XY', URBAN_RESTART_OUT_DTYPE, nohalo=.true.                      ) ! [IN]

       call FILEIO_write( URBAN_SFLX_MW(:,:), basename, URBAN_RESTART_OUT_TITLE,               & ! [IN]
                          VAR_NAME(I_SFLX_MW), VAR_DESC(I_SFLX_MW), VAR_UNIT(I_SFLX_MW),       & ! [IN]
                          'XY', URBAN_RESTART_OUT_DTYPE, nohalo=.true.                         ) ! [IN]
       call FILEIO_write( URBAN_SFLX_MU(:,:), basename, URBAN_RESTART_OUT_TITLE,               & ! [IN]
                          VAR_NAME(I_SFLX_MU), VAR_DESC(I_SFLX_MU), VAR_UNIT(I_SFLX_MU),       & ! [IN]
                          'XY', URBAN_RESTART_OUT_DTYPE, nohalo=.true.                         ) ! [IN]
       call FILEIO_write( URBAN_SFLX_MV(:,:), basename, URBAN_RESTART_OUT_TITLE,               & ! [IN]
                          VAR_NAME(I_SFLX_MV), VAR_DESC(I_SFLX_MV), VAR_UNIT(I_SFLX_MV),       & ! [IN]
                          'XY', URBAN_RESTART_OUT_DTYPE, nohalo=.true.                         ) ! [IN]
       call FILEIO_write( URBAN_SFLX_SH(:,:), basename, URBAN_RESTART_OUT_TITLE,               & ! [IN]
                          VAR_NAME(I_SFLX_SH), VAR_DESC(I_SFLX_SH), VAR_UNIT(I_SFLX_SH),       & ! [IN]
                          'XY', URBAN_RESTART_OUT_DTYPE, nohalo=.true.                         ) ! [IN]
       call FILEIO_write( URBAN_SFLX_LH(:,:), basename, URBAN_RESTART_OUT_TITLE,               & ! [IN]
                          VAR_NAME(I_SFLX_LH), VAR_DESC(I_SFLX_LH), VAR_UNIT(I_SFLX_LH),       & ! [IN]
                          'XY', URBAN_RESTART_OUT_DTYPE, nohalo=.true.                         ) ! [IN]
       call FILEIO_write( URBAN_SFLX_GH(:,:), basename, URBAN_RESTART_OUT_TITLE,               & ! [IN]
                          VAR_NAME(I_SFLX_GH), VAR_DESC(I_SFLX_GH), VAR_UNIT(I_SFLX_GH),       & ! [IN]
                          'XY', URBAN_RESTART_OUT_DTYPE, nohalo=.true.                         ) ! [IN]
       call FILEIO_write( URBAN_SFLX_evap(:,:), basename, URBAN_RESTART_OUT_TITLE,             & ! [IN]
                          VAR_NAME(I_SFLX_evap), VAR_DESC(I_SFLX_evap), VAR_UNIT(I_SFLX_evap), & ! [IN]
                          'XY', URBAN_RESTART_OUT_DTYPE, nohalo=.true.                         ) ! [IN]

    endif

    return
  end subroutine URBAN_vars_restart_write

  !-----------------------------------------------------------------------------
  !> History output set for urban variables
  subroutine URBAN_vars_history
    use scale_history, only: &
       HIST_in
    implicit none
    !---------------------------------------------------------------------------

    if ( URBAN_VARS_CHECKRANGE ) then
       call VALCHECK( URBAN_TR        (IS:IE,JS:JE),          0.0_RP, 1000.0_RP, VAR_NAME(I_TR),        &
                     __FILE__, __LINE__ )
       call VALCHECK( URBAN_TB        (IS:IE,JS:JE),          0.0_RP, 1000.0_RP, VAR_NAME(I_TB),        &
                     __FILE__, __LINE__ )
       call VALCHECK( URBAN_TG        (IS:IE,JS:JE),          0.0_RP, 1000.0_RP, VAR_NAME(I_TG),        &
                     __FILE__, __LINE__ )
       call VALCHECK( URBAN_TC        (IS:IE,JS:JE),          0.0_RP, 1000.0_RP, VAR_NAME(I_TC),        &
                     __FILE__, __LINE__ )
       call VALCHECK( URBAN_QC        (IS:IE,JS:JE),          0.0_RP, 1000.0_RP, VAR_NAME(I_QC),        &
                     __FILE__, __LINE__ )
       call VALCHECK( URBAN_UC        (IS:IE,JS:JE),          0.0_RP, 1000.0_RP, VAR_NAME(I_UC),        &
                     __FILE__, __LINE__ )
       call VALCHECK( URBAN_TRL       (:,IS:IE,JS:JE),        0.0_RP, 1000.0_RP, VAR_NAME(I_TRL),       &
                     __FILE__, __LINE__ )
       call VALCHECK( URBAN_TBL       (:,IS:IE,JS:JE),        0.0_RP, 1000.0_RP, VAR_NAME(I_TBL),       &
                     __FILE__, __LINE__ )
       call VALCHECK( URBAN_TGL       (:,IS:IE,JS:JE),        0.0_RP, 1000.0_RP, VAR_NAME(I_TGL),       &
                     __FILE__, __LINE__ )
       call VALCHECK( URBAN_RAINR     (IS:IE,JS:JE),      -1000.0_RP, 1000.0_RP, VAR_NAME(I_RAINR),     &
                     __FILE__, __LINE__ )
       call VALCHECK( URBAN_RAINB     (IS:IE,JS:JE),      -1000.0_RP, 1000.0_RP, VAR_NAME(I_RAINB),     &
                     __FILE__, __LINE__ )
       call VALCHECK( URBAN_RAING     (IS:IE,JS:JE),      -1000.0_RP, 1000.0_RP, VAR_NAME(I_RAING),     &
                     __FILE__, __LINE__ )
       call VALCHECK( URBAN_ROFF      (IS:IE,JS:JE),      -1000.0_RP, 1000.0_RP, VAR_NAME(I_ROFF),      &
                     __FILE__, __LINE__ )
       call VALCHECK( URBAN_SFC_TEMP  (IS:IE,JS:JE),          0.0_RP, 1000.0_RP, VAR_NAME(I_SFC_TEMP),  &
                     __FILE__, __LINE__ )
       call VALCHECK( URBAN_SFC_albedo(IS:IE,JS:JE,I_LW),     0.0_RP,    2.0_RP, VAR_NAME(I_ALB_LW),    &
                     __FILE__, __LINE__ )
       call VALCHECK( URBAN_SFC_albedo(IS:IE,JS:JE,I_SW),     0.0_RP,    2.0_RP, VAR_NAME(I_ALB_SW),    &
                     __FILE__, __LINE__ )
    endif

    call HIST_in( URBAN_TR(:,:), VAR_NAME(I_TR), VAR_DESC(I_TR), VAR_UNIT(I_TR) )
    call HIST_in( URBAN_TB(:,:), VAR_NAME(I_TB), VAR_DESC(I_TB), VAR_UNIT(I_TB) )
    call HIST_in( URBAN_TG(:,:), VAR_NAME(I_TG), VAR_DESC(I_TG), VAR_UNIT(I_TG) )
    call HIST_in( URBAN_TC(:,:), VAR_NAME(I_TC), VAR_DESC(I_TC), VAR_UNIT(I_TC) )
    call HIST_in( URBAN_QC(:,:), VAR_NAME(I_QC), VAR_DESC(I_QC), VAR_UNIT(I_QC) )
    call HIST_in( URBAN_UC(:,:), VAR_NAME(I_UC), VAR_DESC(I_UC), VAR_UNIT(I_UC) )

    call HIST_in( URBAN_TRL(:,:,:), VAR_NAME(I_TRL), VAR_DESC(I_TRL), VAR_UNIT(I_TRL), zdim='urban' )
    call HIST_in( URBAN_TBL(:,:,:), VAR_NAME(I_TBL), VAR_DESC(I_TBL), VAR_UNIT(I_TBL), zdim='urban' )
    call HIST_in( URBAN_TGL(:,:,:), VAR_NAME(I_TGL), VAR_DESC(I_TGL), VAR_UNIT(I_TGL), zdim='urban' )

    call HIST_in( URBAN_RAINR(:,:), VAR_NAME(I_RAINR), VAR_DESC(I_RAINR), VAR_UNIT(I_RAINR) )
    call HIST_in( URBAN_RAINB(:,:), VAR_NAME(I_RAINB), VAR_DESC(I_RAINB), VAR_UNIT(I_RAINB) )
    call HIST_in( URBAN_RAING(:,:), VAR_NAME(I_RAING), VAR_DESC(I_RAING), VAR_UNIT(I_RAING) )
    call HIST_in( URBAN_ROFF (:,:), VAR_NAME(I_ROFF),  VAR_DESC(I_ROFF),  VAR_UNIT(I_ROFF)  )

    call HIST_in( URBAN_SFC_TEMP  (:,:),      VAR_NAME(I_SFC_TEMP), VAR_DESC(I_SFC_TEMP), VAR_UNIT(I_SFC_TEMP) )
    call HIST_in( URBAN_SFC_albedo(:,:,I_LW), VAR_NAME(I_ALB_LW),   VAR_DESC(I_ALB_LW),   VAR_UNIT(I_ALB_LW)   )
    call HIST_in( URBAN_SFC_albedo(:,:,I_SW), VAR_NAME(I_ALB_SW),   VAR_DESC(I_ALB_SW),   VAR_UNIT(I_ALB_SW)   )

    call HIST_in( URBAN_SFLX_MW  (:,:), VAR_NAME(I_SFLX_MW),   VAR_DESC(I_SFLX_MW),   VAR_UNIT(I_SFLX_MW)   )
    call HIST_in( URBAN_SFLX_MU  (:,:), VAR_NAME(I_SFLX_MU),   VAR_DESC(I_SFLX_MU),   VAR_UNIT(I_SFLX_MU)   )
    call HIST_in( URBAN_SFLX_MV  (:,:), VAR_NAME(I_SFLX_MV),   VAR_DESC(I_SFLX_MV),   VAR_UNIT(I_SFLX_MV)   )
    call HIST_in( URBAN_SFLX_SH  (:,:), VAR_NAME(I_SFLX_SH),   VAR_DESC(I_SFLX_SH),   VAR_UNIT(I_SFLX_SH)   )
    call HIST_in( URBAN_SFLX_LH  (:,:), VAR_NAME(I_SFLX_LH),   VAR_DESC(I_SFLX_LH),   VAR_UNIT(I_SFLX_LH)   )
    call HIST_in( URBAN_SFLX_GH  (:,:), VAR_NAME(I_SFLX_GH),   VAR_DESC(I_SFLX_GH),   VAR_UNIT(I_SFLX_GH)   )
    call HIST_in( URBAN_SFLX_evap(:,:), VAR_NAME(I_SFLX_evap), VAR_DESC(I_SFLX_evap), VAR_UNIT(I_SFLX_evap) )

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
       call STAT_total( total, URBAN_TR(:,:), VAR_NAME(I_TR) )
       call STAT_total( total, URBAN_TB(:,:), VAR_NAME(I_TB) )
       call STAT_total( total, URBAN_TG(:,:), VAR_NAME(I_TG) )
       call STAT_total( total, URBAN_TC(:,:), VAR_NAME(I_TC) )
       call STAT_total( total, URBAN_QC(:,:), VAR_NAME(I_QC) )
       call STAT_total( total, URBAN_UC(:,:), VAR_NAME(I_UC) )

       do k = UKS, UKE
          call STAT_total( total, URBAN_TRL(k,:,:), VAR_NAME(I_TRL) )
          call STAT_total( total, URBAN_TBL(k,:,:), VAR_NAME(I_TBL) )
          call STAT_total( total, URBAN_TGL(k,:,:), VAR_NAME(I_TGL) )
       enddo

       call STAT_total( total, URBAN_RAINR(:,:), VAR_NAME(I_RAINR) )
       call STAT_total( total, URBAN_RAINB(:,:), VAR_NAME(I_RAINB) )
       call STAT_total( total, URBAN_RAING(:,:), VAR_NAME(I_RAING) )
       call STAT_total( total, URBAN_ROFF (:,:), VAR_NAME(I_ROFF)  )

       call STAT_total( total, URBAN_SFC_TEMP  (:,:),      VAR_NAME(I_SFC_TEMP) )
       call STAT_total( total, URBAN_SFC_albedo(:,:,I_LW), VAR_NAME(I_ALB_LW)   )
       call STAT_total( total, URBAN_SFC_albedo(:,:,I_SW), VAR_NAME(I_ALB_SW)   )

       call STAT_total( total, URBAN_SFLX_MW  (:,:), VAR_NAME(I_SFLX_MW)   )
       call STAT_total( total, URBAN_SFLX_MU  (:,:), VAR_NAME(I_SFLX_MU)   )
       call STAT_total( total, URBAN_SFLX_MV  (:,:), VAR_NAME(I_SFLX_MV)   )
       call STAT_total( total, URBAN_SFLX_SH  (:,:), VAR_NAME(I_SFLX_SH)   )
       call STAT_total( total, URBAN_SFLX_LH  (:,:), VAR_NAME(I_SFLX_LH)   )
       call STAT_total( total, URBAN_SFLX_GH  (:,:), VAR_NAME(I_SFLX_GH)   )
       call STAT_total( total, URBAN_SFLX_evap(:,:), VAR_NAME(I_SFLX_evap) )

    endif

    return
  end subroutine URBAN_vars_total

  !-----------------------------------------------------------------------------
  !> Input from External I/O
  subroutine URBAN_vars_external_in( &
       URBAN_TC_in,        &
       URBAN_QC_in,        &
       URBAN_UC_in,        &
       URBAN_SFC_TEMP_in,  &
       URBAN_SFC_albedo_in )
    implicit none

    real(RP), intent(in) :: URBAN_TC_in        (IA,JA)
    real(RP), intent(in) :: URBAN_QC_in        (IA,JA)
    real(RP), intent(in) :: URBAN_UC_in        (IA,JA)
    real(RP), intent(in) :: URBAN_SFC_TEMP_in  (IA,JA)
    real(RP), intent(in) :: URBAN_SFC_albedo_in(IA,JA,2)

    integer :: k
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** External Input (coupler) ***'

    URBAN_TR(:,:) = URBAN_SFC_TEMP_in(:,:)
    URBAN_TB(:,:) = URBAN_SFC_TEMP_in(:,:)
    URBAN_TG(:,:) = URBAN_SFC_TEMP_in(:,:)

    URBAN_TC(:,:) = URBAN_TC_in(:,:)
    URBAN_QC(:,:) = URBAN_QC_in(:,:)
    URBAN_UC(:,:) = URBAN_UC_in(:,:)

    do k = UKS, UKE
       URBAN_TRL(k,:,:) = URBAN_SFC_TEMP_in(:,:)
       URBAN_TBL(k,:,:) = URBAN_SFC_TEMP_in(:,:)
       URBAN_TGL(k,:,:) = URBAN_SFC_TEMP_in(:,:)
    end do

    URBAN_RAINR(:,:) = 0.0_RP
    URBAN_RAINB(:,:) = 0.0_RP
    URBAN_RAING(:,:) = 0.0_RP
    URBAN_ROFF (:,:) = 0.0_RP

    URBAN_SFC_TEMP  (:,:)   = URBAN_SFC_TEMP_in  (:,:)
    URBAN_SFC_albedo(:,:,:) = URBAN_SFC_albedo_in(:,:,:)

    URBAN_SFLX_MW  (:,:) = 0.0_RP
    URBAN_SFLX_MU  (:,:) = 0.0_RP
    URBAN_SFLX_MV  (:,:) = 0.0_RP
    URBAN_SFLX_SH  (:,:) = 0.0_RP
    URBAN_SFLX_LH  (:,:) = 0.0_RP
    URBAN_SFLX_GH  (:,:) = 0.0_RP
    URBAN_SFLX_evap(:,:) = 0.0_RP

    call URBAN_vars_total

    return
  end subroutine URBAN_vars_external_in

end module mod_urban_vars
