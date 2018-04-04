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
!! @li      2014-08-28 (A.Noda)     [mod] input of ATMOS_PHY_SF_SFC_{TEMP,albedo}
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
  use scale_atmos_grid_cartesC_index
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
  public :: ATMOS_PHY_SF_vars_external_in

  public :: ATMOS_PHY_SF_vars_restart_create
  public :: ATMOS_PHY_SF_vars_restart_open
  public :: ATMOS_PHY_SF_vars_restart_def_var
  public :: ATMOS_PHY_SF_vars_restart_enddef
  public :: ATMOS_PHY_SF_vars_restart_close

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical,               public :: ATMOS_PHY_SF_RESTART_OUTPUT                 = .false.                !< output restart file?

  character(len=H_LONG),  public :: ATMOS_PHY_SF_RESTART_IN_BASENAME           = ''                     !< Basename of the input  file
  logical,                public :: ATMOS_PHY_SF_RESTART_IN_AGGREGATE                                   !< Switch to use aggregate file
  logical,                public :: ATMOS_PHY_SF_RESTART_IN_POSTFIX_TIMELABEL  = .false.                !< Add timelabel to the basename of input  file?
  character(len=H_LONG),  public :: ATMOS_PHY_SF_RESTART_OUT_BASENAME          = ''                     !< Basename of the output file
  logical,                public :: ATMOS_PHY_SF_RESTART_OUT_AGGREGATE                                  !< Switch to use aggregate file
  logical,                public :: ATMOS_PHY_SF_RESTART_OUT_POSTFIX_TIMELABEL = .true.                 !< Add timelabel to the basename of output file?
  character(len=H_MID),   public :: ATMOS_PHY_SF_RESTART_OUT_TITLE             = 'ATMOS_PHY_SF restart' !< title    of the output file
  character(len=H_SHORT), public :: ATMOS_PHY_SF_RESTART_OUT_DTYPE             = 'DEFAULT'              !< REAL4 or REAL8

  real(RP), public, allocatable :: ATMOS_PHY_SF_DENS_t(:,:)   ! tendency DENS [    kg/m3/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_MOMZ_t(:,:)   ! tendency MOMZ [m/s*kg/m3/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_RHOU_t(:,:)   ! tendency rho*U [m/s*kg/m3/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_RHOV_t(:,:)   ! tendency rho*V [m/s*kg/m3/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_RHOH  (:,:)   ! diabatic heating [J/m3/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_RHOT_t(:,:)   ! tendency RHOT [K  *kg/m3/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_RHOQ_t(:,:,:) ! tendency rho*QTRC [    kg/kg/s]

  real(RP), public, allocatable :: ATMOS_PHY_SF_SFC_TEMP  (:,:)   ! surface skin temperature             [K]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFC_albedo(:,:,:) ! surface albedo                       (0-1)
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFC_Z0M   (:,:)   ! surface roughness length, ocean only [m]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFC_Z0H   (:,:)   ! surface roughness length, ocean only [m]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFC_Z0E   (:,:)   ! surface roughness length, ocean only [m]

  ! surface diagnostic variables
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFC_DENS  (:,:)   ! surface atmosphere density  [kg/m3]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFC_PRES  (:,:)   ! surface atmosphere pressure [Pa]

  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_MW   (:,:)   ! z-momentum flux (area center) [m/s*kg/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_MU   (:,:)   ! x-momentum flux (area center) [m/s*kg/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_MV   (:,:)   ! y-momentum flux (area center) [m/s*kg/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_SH   (:,:)   ! sensible heat flux [J/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_LH   (:,:)   ! latent heat flux [J/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_GH   (:,:)   ! ground heat flux [J/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_SFLX_QTRC (:,:,:) ! tracer mass flux [kg/m2/s]

  real(RP), public, allocatable :: ATMOS_PHY_SF_U10       (:,:)   ! 10m x-wind [m/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_V10       (:,:)   ! 10m y-wind [m/s]
  real(RP), public, allocatable :: ATMOS_PHY_SF_T2        (:,:)   ! 2m temperature [K]
  real(RP), public, allocatable :: ATMOS_PHY_SF_Q2        (:,:)   ! 2m specific humidity [kg/kg]

  real(RP), public, allocatable :: ATMOS_PHY_SF_l_mo      (:,:)   ! monin-obukov length

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
  integer,                private, parameter :: VMAX            = 6 !< number of the variables
  integer,                private, parameter :: I_SFC_TEMP      = 1
  integer,                private, parameter :: I_SFC_albedo_LW = 2
  integer,                private, parameter :: I_SFC_albedo_SW = 3
  integer,                private, parameter :: I_SFC_Z0M       = 4
  integer,                private, parameter :: I_SFC_Z0H       = 5
  integer,                private, parameter :: I_SFC_Z0E       = 6

  character(len=H_SHORT), private            :: VAR_NAME(VMAX) !< name  of the variables
  character(len=H_MID),   private            :: VAR_DESC(VMAX) !< desc. of the variables
  character(len=H_MID),   private            :: VAR_STDN(VMAX) !< standard name of the variables
  character(len=H_SHORT), private            :: VAR_UNIT(VMAX) !< unit  of the variables
  integer,                private            :: VAR_ID(VMAX)   !< ID    of the variables
  integer,                private            :: restart_fid = -1  ! file ID

  data VAR_NAME / 'SFC_TEMP',   &
                  'SFC_ALB_LW', &
                  'SFC_ALB_SW', &
                  'SFC_Z0M',    &
                  'SFC_Z0H',    &
                  'SFC_Z0E'     /

  data VAR_DESC / 'surface skin temperature',            &
                  'surface albedo for longwave',         &
                  'surface albedo for shortwave',        &
                  'surface roughness length (momentum)', &
                  'surface roughness length (heat)',     &
                  'surface roughness length (vapor)'     /

  data VAR_STDN / 'surface_temp', &
                  '', &
                  '', &
                  'surface_roughness_length_for_momentum_in_air', &
                  'surface_roughness_length_for_heat_in_air', &
                  '' /

  data VAR_UNIT / 'K', &
                  '1', &
                  '1', &
                  'm', &
                  'm', &
                  'm'  /

  real(RP), private :: ATMOS_PHY_SF_DEFAULT_SFC_TEMP   = 300.0_RP
  real(RP), private :: ATMOS_PHY_SF_DEFAULT_SFC_albedo = 0.4_RP
  real(RP), private :: ATMOS_PHY_SF_DEFAULT_SFC_Z0     = 1E-4_RP

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_SF_vars_setup
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    implicit none

    NAMELIST / PARAM_ATMOS_PHY_SF_VARS / &
       ATMOS_PHY_SF_RESTART_IN_BASENAME,           &
       ATMOS_PHY_SF_RESTART_IN_AGGREGATE,          &
       ATMOS_PHY_SF_RESTART_IN_POSTFIX_TIMELABEL,  &
       ATMOS_PHY_SF_RESTART_OUTPUT,                &
       ATMOS_PHY_SF_RESTART_OUT_BASENAME,          &
       ATMOS_PHY_SF_RESTART_OUT_AGGREGATE,         &
       ATMOS_PHY_SF_RESTART_OUT_POSTFIX_TIMELABEL, &
       ATMOS_PHY_SF_RESTART_OUT_TITLE,             &
       ATMOS_PHY_SF_RESTART_OUT_DTYPE,             &
       ATMOS_PHY_SF_DEFAULT_SFC_TEMP,              &
       ATMOS_PHY_SF_DEFAULT_SFC_albedo,            &
       ATMOS_PHY_SF_DEFAULT_SFC_Z0

    integer :: ierr
    integer :: iv
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[VARS] / Categ[ATMOS PHY_SF] / Origin[SCALE-RM]'

    allocate( ATMOS_PHY_SF_DENS_t(IA,JA)    )
    allocate( ATMOS_PHY_SF_MOMZ_t(IA,JA)    )
    allocate( ATMOS_PHY_SF_RHOU_t(IA,JA)    )
    allocate( ATMOS_PHY_SF_RHOV_t(IA,JA)    )
    allocate( ATMOS_PHY_SF_RHOT_t(IA,JA)    )
    allocate( ATMOS_PHY_SF_RHOH  (IA,JA)    )
    allocate( ATMOS_PHY_SF_RHOQ_t(IA,JA,QA) )
    ATMOS_PHY_SF_DENS_t(:,:)   = UNDEF
    ATMOS_PHY_SF_MOMZ_t(:,:)   = UNDEF
    ATMOS_PHY_SF_RHOU_t(:,:)   = UNDEF
    ATMOS_PHY_SF_RHOV_t(:,:)   = UNDEF
    ATMOS_PHY_SF_RHOT_t(:,:)   = UNDEF
    ATMOS_PHY_SF_RHOH  (:,:)   = UNDEF
    ATMOS_PHY_SF_RHOQ_t(:,:,:) = UNDEF

    allocate( ATMOS_PHY_SF_SFC_TEMP  (IA,JA)    )
    allocate( ATMOS_PHY_SF_SFC_albedo(IA,JA,2)  )
    allocate( ATMOS_PHY_SF_SFC_Z0M   (IA,JA)    )
    allocate( ATMOS_PHY_SF_SFC_Z0H   (IA,JA)    )
    allocate( ATMOS_PHY_SF_SFC_Z0E   (IA,JA)    )
    ATMOS_PHY_SF_SFC_TEMP  (:,:)   = UNDEF
    ATMOS_PHY_SF_SFC_albedo(:,:,:) = UNDEF
    ATMOS_PHY_SF_SFC_Z0M   (:,:)   = UNDEF
    ATMOS_PHY_SF_SFC_Z0H   (:,:)   = UNDEF
    ATMOS_PHY_SF_SFC_Z0E   (:,:)   = UNDEF

    allocate( ATMOS_PHY_SF_SFC_DENS  (IA,JA)    )
    allocate( ATMOS_PHY_SF_SFC_PRES  (IA,JA)    )
    ATMOS_PHY_SF_SFC_DENS  (:,:)   = UNDEF
    ATMOS_PHY_SF_SFC_PRES  (:,:)   = UNDEF

    allocate( ATMOS_PHY_SF_SFLX_MW   (IA,JA)    )
    allocate( ATMOS_PHY_SF_SFLX_MU   (IA,JA)    )
    allocate( ATMOS_PHY_SF_SFLX_MV   (IA,JA)    )
    allocate( ATMOS_PHY_SF_SFLX_SH   (IA,JA)    )
    allocate( ATMOS_PHY_SF_SFLX_LH   (IA,JA)    )
    allocate( ATMOS_PHY_SF_SFLX_GH   (IA,JA)    )
    allocate( ATMOS_PHY_SF_SFLX_QTRC (IA,JA,max(QA,1)) )
    ATMOS_PHY_SF_SFLX_MW  (:,:)   = UNDEF
    ATMOS_PHY_SF_SFLX_MU  (:,:)   = UNDEF
    ATMOS_PHY_SF_SFLX_MV  (:,:)   = UNDEF
    ATMOS_PHY_SF_SFLX_SH  (:,:)   = UNDEF
    ATMOS_PHY_SF_SFLX_LH  (:,:)   = UNDEF
    ATMOS_PHY_SF_SFLX_GH  (:,:)   = UNDEF
    ATMOS_PHY_SF_SFLX_QTRC(:,:,:) = UNDEF

    allocate( ATMOS_PHY_SF_U10       (IA,JA)    )
    allocate( ATMOS_PHY_SF_V10       (IA,JA)    )
    allocate( ATMOS_PHY_SF_T2        (IA,JA)    )
    allocate( ATMOS_PHY_SF_Q2        (IA,JA)    )
    allocate( ATMOS_PHY_SF_l_mo      (IA,JA)    )
    ATMOS_PHY_SF_U10      (:,:)   = UNDEF
    ATMOS_PHY_SF_V10      (:,:)   = UNDEF
    ATMOS_PHY_SF_T2       (:,:)   = UNDEF
    ATMOS_PHY_SF_Q2       (:,:)   = UNDEF
    ATMOS_PHY_SF_l_mo     (:,:)   = UNDEF

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_SF_VARS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_SF_VARS. Check!'
       call PRC_abort
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_ATMOS_PHY_SF_VARS)

    ! [add] 2014/08/28 A.Noda
    ATMOS_PHY_SF_SFC_TEMP  (:,:)   = ATMOS_PHY_SF_DEFAULT_SFC_TEMP
    ATMOS_PHY_SF_SFC_albedo(:,:,:) = ATMOS_PHY_SF_DEFAULT_SFC_albedo
    ATMOS_PHY_SF_SFC_Z0M   (:,:)   = ATMOS_PHY_SF_DEFAULT_SFC_Z0
    ATMOS_PHY_SF_SFC_Z0H   (:,:)   = ATMOS_PHY_SF_DEFAULT_SFC_Z0
    ATMOS_PHY_SF_SFC_Z0E   (:,:)   = ATMOS_PHY_SF_DEFAULT_SFC_Z0

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** [ATMOS_PHY_SF] prognostic/diagnostic variables'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A24,A,A48,A,A12,A)') &
               '***       |', 'VARNAME                 ','|', &
               'DESCRIPTION                                     ', '[', 'UNIT        ', ']'
    do iv = 1, VMAX
       if( IO_L ) write(IO_FID_LOG,'(1x,A,I3,A,A24,A,A48,A,A12,A)') &
                  '*** NO.',iv,'|',VAR_NAME(iv),'|',VAR_DESC(iv),'[',VAR_UNIT(iv),']'
    enddo

    if( IO_L ) write(IO_FID_LOG,*)
    if ( ATMOS_PHY_SF_RESTART_IN_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart input?  : YES, file = ', trim(ATMOS_PHY_SF_RESTART_IN_BASENAME)
       if( IO_L ) write(IO_FID_LOG,*) '*** Add timelabel?  : ', ATMOS_PHY_SF_RESTART_IN_POSTFIX_TIMELABEL
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart input?  : NO'
    endif
    if (       ATMOS_PHY_SF_RESTART_OUTPUT             &
         .AND. ATMOS_PHY_SF_RESTART_OUT_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart output? : YES, file = ', trim(ATMOS_PHY_SF_RESTART_OUT_BASENAME)
       if( IO_L ) write(IO_FID_LOG,*) '*** Add timelabel?  : ', ATMOS_PHY_SF_RESTART_OUT_POSTFIX_TIMELABEL
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

    call COMM_vars8( ATMOS_PHY_SF_SFC_TEMP  (:,:),       1 )
    call COMM_vars8( ATMOS_PHY_SF_SFC_albedo(:,:,I_LW),  2 )
    call COMM_vars8( ATMOS_PHY_SF_SFC_albedo(:,:,I_SW),  3 )
    call COMM_vars8( ATMOS_PHY_SF_SFC_Z0M   (:,:),       4 )
    call COMM_vars8( ATMOS_PHY_SF_SFC_Z0H   (:,:),       5 )
    call COMM_vars8( ATMOS_PHY_SF_SFC_Z0E   (:,:),       6 )

    call COMM_wait ( ATMOS_PHY_SF_SFC_TEMP  (:,:),       1 )
    call COMM_wait ( ATMOS_PHY_SF_SFC_albedo(:,:,I_LW),  2 )
    call COMM_wait ( ATMOS_PHY_SF_SFC_albedo(:,:,I_SW),  3 )
    call COMM_wait ( ATMOS_PHY_SF_SFC_Z0M   (:,:),       4 )
    call COMM_wait ( ATMOS_PHY_SF_SFC_Z0H   (:,:),       5 )
    call COMM_wait ( ATMOS_PHY_SF_SFC_Z0E   (:,:),       6 )

    return
  end subroutine ATMOS_PHY_SF_vars_fillhalo

  !-----------------------------------------------------------------------------
  !> Open restart file for read
  subroutine ATMOS_PHY_SF_vars_restart_open
    use scale_time, only: &
       TIME_gettimelabel
    use scale_file_cartesC, only: &
       FILE_CARTESC_open
    implicit none

    character(len=19)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Open restart file (ATMOS_PHY_SF) ***'

    if ( ATMOS_PHY_SF_RESTART_IN_BASENAME /= '' ) then

       if ( ATMOS_PHY_SF_RESTART_IN_POSTFIX_TIMELABEL ) then
          call TIME_gettimelabel( timelabel )
          basename = trim(ATMOS_PHY_SF_RESTART_IN_BASENAME)//'_'//trim(timelabel)
       else
          basename = trim(ATMOS_PHY_SF_RESTART_IN_BASENAME)
       endif

       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(basename)

       call FILE_CARTESC_open( basename, restart_fid, aggregate=ATMOS_PHY_SF_RESTART_IN_AGGREGATE )

    else
       if( IO_L ) write(IO_FID_LOG,*) '*** restart file for ATMOS_PHY_SF is not specified.'
       ATMOS_PHY_SF_SFC_TEMP  (:,:)   = ATMOS_PHY_SF_DEFAULT_SFC_TEMP
       ATMOS_PHY_SF_SFC_albedo(:,:,:) = ATMOS_PHY_SF_DEFAULT_SFC_albedo
       ATMOS_PHY_SF_SFC_Z0M   (:,:)   = ATMOS_PHY_SF_DEFAULT_SFC_Z0
       ATMOS_PHY_SF_SFC_Z0H   (:,:)   = ATMOS_PHY_SF_DEFAULT_SFC_Z0
       ATMOS_PHY_SF_SFC_Z0E   (:,:)   = ATMOS_PHY_SF_DEFAULT_SFC_Z0
    endif

    return
  end subroutine ATMOS_PHY_SF_vars_restart_open

  !-----------------------------------------------------------------------------
  !> Read restart
  subroutine ATMOS_PHY_SF_vars_restart_read
    use scale_const, only: &
       I_SW => CONST_I_SW, &
       I_LW => CONST_I_LW
    use scale_file, only: &
       FILE_get_aggregate
    use scale_file_cartesC, only: &
       FILE_CARTESC_read, &
       FILE_CARTESC_flush
    implicit none
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Read from restart file (ATMOS_PHY_SF) ***'

       call FILE_CARTESC_read( restart_fid, VAR_NAME(1), 'XY', & ! [IN]
                               ATMOS_PHY_SF_SFC_TEMP  (:,:)    ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(2), 'XY',   & ! [IN]
                               ATMOS_PHY_SF_SFC_albedo(:,:,I_LW) ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(3), 'XY',   & ! [IN]
                               ATMOS_PHY_SF_SFC_albedo(:,:,I_SW) ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(4), 'XY', & ! [IN]
                               ATMOS_PHY_SF_SFC_Z0M   (:,:)    ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(5), 'XY', & ! [IN]
                               ATMOS_PHY_SF_SFC_Z0H   (:,:)    ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(6), 'XY', & ! [IN]
                               ATMOS_PHY_SF_SFC_Z0E   (:,:)    ) ! [OUT]

       if ( FILE_get_AGGREGATE(restart_fid) ) then
          call FILE_CARTESC_flush( restart_fid ) ! X/Y halos have been read from file
       else
          call ATMOS_PHY_SF_vars_fillhalo
       end if

       call ATMOS_PHY_SF_vars_checktotal

    else
       if( IO_L ) write(IO_FID_LOG,*) '*** invalid restart file ID for ATMOS_PHY_SF.'
    endif

    return
  end subroutine ATMOS_PHY_SF_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Input from External I/O
  subroutine ATMOS_PHY_SF_vars_external_in( &
      sfc_temp_in,   &
      albedo_lw_in,  &
      albedo_sw_in,  &
      sfc_z0m_in,    &
      sfc_z0h_in,    &
      sfc_z0e_in     )
    use scale_const, only: &
       I_SW => CONST_I_SW, &
       I_LW => CONST_I_LW
    implicit none

    real(RP), intent(in) :: sfc_temp_in (:,:)
    real(RP), intent(in) :: albedo_lw_in(:,:)
    real(RP), intent(in) :: albedo_sw_in(:,:)
    real(RP), intent(in) :: sfc_z0m_in  (:,:)
    real(RP), intent(in) :: sfc_z0h_in  (:,:)
    real(RP), intent(in) :: sfc_z0e_in  (:,:)
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** External Input (PHY_SF) ***'

    ATMOS_PHY_SF_SFC_TEMP   (:,:)      = sfc_temp_in (:,:)
    ATMOS_PHY_SF_SFC_albedo (:,:,I_LW) = albedo_lw_in(:,:)
    ATMOS_PHY_SF_SFC_albedo (:,:,I_SW) = albedo_sw_in(:,:)
    ATMOS_PHY_SF_SFC_Z0M    (:,:)      = sfc_z0m_in  (:,:)
    ATMOS_PHY_SF_SFC_Z0H    (:,:)      = sfc_z0h_in  (:,:)
    ATMOS_PHY_SF_SFC_Z0E    (:,:)      = sfc_z0e_in  (:,:)

    call ATMOS_PHY_SF_vars_fillhalo

    return
  end subroutine ATMOS_PHY_SF_vars_external_in

  !-----------------------------------------------------------------------------
  !> Create restart file
  subroutine ATMOS_PHY_SF_vars_restart_create
    use scale_time, only: &
       TIME_gettimelabel
    use scale_file_cartesC, only: &
       FILE_CARTESC_create
    implicit none

    character(len=19)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    if ( ATMOS_PHY_SF_RESTART_OUT_BASENAME /= '' ) then

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Create restart file (ATMOS_PHY_AE) ***'

       if ( ATMOS_PHY_SF_RESTART_OUT_POSTFIX_TIMELABEL ) then
          call TIME_gettimelabel( timelabel )
          basename = trim(ATMOS_PHY_SF_RESTART_OUT_BASENAME)//'_'//trim(timelabel)
       else
          basename = trim(ATMOS_PHY_SF_RESTART_OUT_BASENAME)
       endif

       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(basename)

       call FILE_CARTESC_create( &
            basename, ATMOS_PHY_SF_RESTART_OUT_TITLE, ATMOS_PHY_SF_RESTART_OUT_DTYPE, & ! [IN]
            restart_fid,                                                              & ! [OUT]
            aggregate=ATMOS_PHY_SF_RESTART_OUT_AGGREGATE                              ) ! [IN]
    endif

    return
  end subroutine ATMOS_PHY_SF_vars_restart_create

  !-----------------------------------------------------------------------------
  !> Exit netCDF define mode
  subroutine ATMOS_PHY_SF_vars_restart_enddef
    use scale_file_cartesC, only: &
       FILE_CARTESC_enddef
    implicit none

    if ( restart_fid /= -1 ) then
       call FILE_CARTESC_enddef( restart_fid ) ! [IN]
    endif

    return
  end subroutine ATMOS_PHY_SF_vars_restart_enddef

  !-----------------------------------------------------------------------------
  !> Close restart file
  subroutine ATMOS_PHY_SF_vars_restart_close
    use scale_file_cartesC, only: &
       FILE_CARTESC_close
    implicit none
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Close restart file (ATMOS_PHY_SF) ***'

       call FILE_CARTESC_close( restart_fid ) ! [IN]

       restart_fid = -1
    endif

    return
  end subroutine ATMOS_PHY_SF_vars_restart_close

  !-----------------------------------------------------------------------------
  !> Write restart
  subroutine ATMOS_PHY_SF_vars_restart_def_var
    use scale_file_cartesC, only: &
       FILE_CARTESC_def_var
    implicit none
    integer :: i
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then

       do i = 1, 6
          call FILE_CARTESC_def_var( restart_fid,     & ! [IN]
               VAR_NAME(i), VAR_DESC(i), VAR_UNIT(i), & ! [IN]
               'XY', ATMOS_PHY_SF_RESTART_OUT_DTYPE,  & ! [IN]
               VAR_ID(i),                             & ! [OUT]
               standard_name=VAR_STDN(i)              ) ! [IN]
       end do

    endif

    return
  end subroutine ATMOS_PHY_SF_vars_restart_def_var

  !-----------------------------------------------------------------------------
  !> Write variables to restart file
  subroutine ATMOS_PHY_SF_vars_restart_write
    use scale_const, only: &
       I_SW => CONST_I_SW, &
       I_LW => CONST_I_LW
    use scale_file_cartesC, only: &
       FILE_CARTESC_write_var
    implicit none

    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then

       call ATMOS_PHY_SF_vars_fillhalo

       call ATMOS_PHY_SF_vars_checktotal

       call FILE_CARTESC_write_var( restart_fid, VAR_ID(1), ATMOS_PHY_SF_SFC_TEMP  (:,:),      & ! [IN]
                              VAR_NAME(1), 'XY'                                          ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(2), ATMOS_PHY_SF_SFC_albedo(:,:,I_LW), & ! [IN]
                              VAR_NAME(2), 'XY'                                          ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(3), ATMOS_PHY_SF_SFC_albedo(:,:,I_SW), & ! [IN]
                              VAR_NAME(3), 'XY'                                          ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(4), ATMOS_PHY_SF_SFC_Z0M   (:,:),      & ! [IN]
                              VAR_NAME(4), 'XY'                                          ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(5), ATMOS_PHY_SF_SFC_Z0H   (:,:),      & ! [IN]
                              VAR_NAME(5), 'XY'                                          ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(6), ATMOS_PHY_SF_SFC_Z0E   (:,:),      & ! [IN]
                              VAR_NAME(6), 'XY'                                          ) ! [IN]

    endif

    return
  end subroutine ATMOS_PHY_SF_vars_restart_write

  subroutine ATMOS_PHY_SF_vars_checktotal
    use scale_const, only: &
       I_SW => CONST_I_SW, &
       I_LW => CONST_I_LW
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STATISTICS_total
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_AREA, &
       ATMOS_GRID_CARTESC_REAL_TOTAREA

    !---------------------------------------------------------------------------

    if ( STATISTICS_checktotal ) then
       call STATISTICS_total( IA, IS, IE, JA, JS, JE, &
                              ATMOS_PHY_SF_SFC_TEMP  (:,:),      VAR_NAME(1), & ! (in)
                              ATMOS_GRID_CARTESC_REAL_AREA(:,:),              & ! (in)
                              ATMOS_GRID_CARTESC_REAL_TOTAREA                 ) ! (in)
       call STATISTICS_total( IA, IS, IE, JA, JS, JE, &
                              ATMOS_PHY_SF_SFC_albedo(:,:,I_LW), VAR_NAME(2), & ! (in)
                              ATMOS_GRID_CARTESC_REAL_AREA(:,:),              & ! (in)
                              ATMOS_GRID_CARTESC_REAL_TOTAREA                 ) ! (in)
       call STATISTICS_total( IA, IS, IE, JA, JS, JE, &
                              ATMOS_PHY_SF_SFC_albedo(:,:,I_SW), VAR_NAME(3), & ! (in)
                              ATMOS_GRID_CARTESC_REAL_AREA(:,:),              & ! (in)
                              ATMOS_GRID_CARTESC_REAL_TOTAREA                 ) ! (in)
       call STATISTICS_total( IA, IS, IE, JA, JS, JE, &
                              ATMOS_PHY_SF_SFC_Z0M   (:,:),      VAR_NAME(4), & ! (in)
                              ATMOS_GRID_CARTESC_REAL_AREA(:,:),              & ! (in)
                              ATMOS_GRID_CARTESC_REAL_TOTAREA                 ) ! (in)
       call STATISTICS_total( IA, IS, IE, JA, JS, JE, &
                              ATMOS_PHY_SF_SFC_Z0H   (:,:),      VAR_NAME(5), & ! (in)
                              ATMOS_GRID_CARTESC_REAL_AREA(:,:),              & ! (in)
                              ATMOS_GRID_CARTESC_REAL_TOTAREA                 ) ! (in)
       call STATISTICS_total( IA, IS, IE, JA, JS, JE, &
                              ATMOS_PHY_SF_SFC_Z0E   (:,:),      VAR_NAME(6), & ! (in)
                              ATMOS_GRID_CARTESC_REAL_AREA(:,:),              & ! (in)
                              ATMOS_GRID_CARTESC_REAL_TOTAREA                 ) ! (in)
    endif

    return
  end subroutine ATMOS_PHY_SF_vars_checktotal

end module mod_atmos_phy_sf_vars
