!-------------------------------------------------------------------------------
!> module Atmosphere / Physics Cumulus
!!
!! @par Description
!!          Container for mod_atmos_phy_cp
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_atmos_phy_cp_vars
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_debug
  use scale_atmos_grid_cartesC_index
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_CP_vars_setup
  public :: ATMOS_PHY_CP_vars_fillhalo
  public :: ATMOS_PHY_CP_vars_restart_read
  public :: ATMOS_PHY_CP_vars_restart_write

  public :: ATMOS_PHY_CP_vars_restart_create
  public :: ATMOS_PHY_CP_vars_restart_open
  public :: ATMOS_PHY_CP_vars_restart_def_var
  public :: ATMOS_PHY_CP_vars_restart_enddef
  public :: ATMOS_PHY_CP_vars_restart_close

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical,               public :: ATMOS_PHY_CP_RESTART_OUTPUT                 = .false.                !< output restart file?

  character(len=H_LONG),  public :: ATMOS_PHY_CP_RESTART_IN_BASENAME           = ''                     !< Basename of the input  file
  logical,                public :: ATMOS_PHY_CP_RESTART_IN_POSTFIX_TIMELABEL  = .false.                !< Add timelabel to the basename of input  file?
  logical,                public :: ATMOS_PHY_CP_RESTART_IN_AGGREGATE                                   !< Switch to use aggregate file
  character(len=H_LONG),  public :: ATMOS_PHY_CP_RESTART_OUT_BASENAME          = ''                     !< Basename of the output file
  logical,                public :: ATMOS_PHY_CP_RESTART_OUT_AGGREGATE                                  !< Switch to use aggregate file
  logical,                public :: ATMOS_PHY_CP_RESTART_OUT_POSTFIX_TIMELABEL = .true.                 !< Add timelabel to the basename of output file?
  character(len=H_MID),   public :: ATMOS_PHY_CP_RESTART_OUT_TITLE             = 'ATMOS_PHY_CP restart' !< title    of the output file
  character(len=H_SHORT), public :: ATMOS_PHY_CP_RESTART_OUT_DTYPE             = 'DEFAULT'              !< REAL4 or REAL8

  ! restart variables
  real(RP), public, allocatable :: ATMOS_PHY_CP_DENS_t  (:,:,:)   ! tendency DENS [kg/m3/s]
  real(RP), public, allocatable :: ATMOS_PHY_CP_MOMZ_t  (:,:,:)   ! tendency MOMZ [kg/m2/s2]
  real(RP), public, allocatable :: ATMOS_PHY_CP_RHOT_t  (:,:,:)   ! tendency RHOT [K*kg/m3/s]
  real(RP), public, allocatable :: ATMOS_PHY_CP_RHOQV_t (:,:,:)   ! tendency rho*QV   [kg/kg/s]
  real(RP), public, allocatable :: ATMOS_PHY_CP_RHOHYD_t(:,:,:,:) ! tendency rho*QHYD [kg/kg/s]

  ! only for K-F scheme
  real(RP), public, allocatable :: ATMOS_PHY_CP_w0mean        (:,:,:) ! running mean vertical wind velocity [m/s]
  real(RP), public, allocatable :: ATMOS_PHY_CP_kf_nca        (:,:)   ! advection/cumulus convection timescale for KF[sec]


  ! diagnostic variables
  real(RP), public, allocatable :: ATMOS_PHY_CP_MFLX_cloudbase(:,:)   ! cloud base mass flux [kg/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_CP_SFLX_rain     (:,:)   ! convective rain [kg/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_CP_SFLX_snow     (:,:)   ! convective snow [kg/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_CP_SFLX_ENGI     (:,:)   ! internal energy [J/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_CP_cloudtop      (:,:)   ! cloud top  height [m]
  real(RP), public, allocatable :: ATMOS_PHY_CP_cloudbase     (:,:)   ! cloud base height [m]
  real(RP), public, allocatable :: ATMOS_PHY_CP_cldfrac_dp    (:,:,:) ! cloud fraction (deep    convection) (0-1)
  real(RP), public, allocatable :: ATMOS_PHY_CP_cldfrac_sh    (:,:,:) ! cloud fraction (shallow convection) (0-1)

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,                private, parameter :: VMAX = 2 !< number of the restart variables
  integer,                private, parameter :: I_w0mean = 1
  integer,                private, parameter :: I_kf_nca = 2


  character(len=H_SHORT), private            :: VAR_NAME(VMAX) !< name  of the variables
  character(len=H_MID),   private            :: VAR_DESC(VMAX) !< desc. of the variables
  character(len=H_SHORT), private            :: VAR_UNIT(VMAX) !< unit  of the variables
  character(len=H_SHORT), private            :: VAR_DIM (VMAX) !< dimension type
  integer,                private            :: VAR_ID  (VMAX) !< ID    of the variables
  integer,                private            :: restart_fid = -1  ! file ID

  data VAR_NAME / 'w0mean', &
                  'kf_nca'  /

  data VAR_DESC / 'running mean vertical velocity',                &
                  'advection/cumulus convection timescale for KF'  /

  data VAR_UNIT / 'm/s', &
                  'sec'  /

  data VAR_DIM  / 'ZXY',  &
                  'XY'    /


  ! tendency names
  integer,                private              :: VMAX_t       !< number of the tendency variables dens+rhoh+qv+N_HYD
  integer,                private              :: I_cp_dens_t = 1
  integer,                private              :: I_cp_rhot_t = 2
  integer,                private              :: I_cp_qv_t   = 3

  character(len=H_SHORT), private, allocatable :: VAR_t_NAME(:) !< name  of the variables
  character(len=H_MID),   private, allocatable :: VAR_t_DESC(:) !< desc. of the variables
  character(len=H_SHORT), private, allocatable :: VAR_t_UNIT(:) !< unit  of the variables
  integer,                private, allocatable :: VAR_t_ID  (:) !< ID    of the variables

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_CP_vars_setup
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_atmos_hydrometeor, only: &
       N_HYD, &
       HYD_NAME
    implicit none

    namelist / PARAM_ATMOS_PHY_CP_VARS / &
       ATMOS_PHY_CP_RESTART_IN_BASENAME,           &
       ATMOS_PHY_CP_RESTART_IN_AGGREGATE,          &
       ATMOS_PHY_CP_RESTART_IN_POSTFIX_TIMELABEL,  &
       ATMOS_PHY_CP_RESTART_OUTPUT,                &
       ATMOS_PHY_CP_RESTART_OUT_BASENAME,          &
       ATMOS_PHY_CP_RESTART_OUT_AGGREGATE,         &
       ATMOS_PHY_CP_RESTART_OUT_POSTFIX_TIMELABEL, &
       ATMOS_PHY_CP_RESTART_OUT_TITLE,             &
       ATMOS_PHY_CP_RESTART_OUT_DTYPE

    integer :: ierr
    integer :: iv
    integer :: iq
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_CP_vars_setup",*) 'Setup'

    allocate( ATMOS_PHY_CP_DENS_t  (KA,IA,JA)       )
    allocate( ATMOS_PHY_CP_MOMZ_t  (KA,IA,JA)       )
    allocate( ATMOS_PHY_CP_RHOT_t  (KA,IA,JA)       )
    allocate( ATMOS_PHY_CP_RHOQV_t (KA,IA,JA)       )
    allocate( ATMOS_PHY_CP_RHOHYD_t(KA,IA,JA,N_HYD) )
    ATMOS_PHY_CP_DENS_t  (:,:,:)   = 0.0_RP
    ATMOS_PHY_CP_MOMZ_t  (:,:,:)   = UNDEF
    ATMOS_PHY_CP_RHOT_t  (:,:,:)   = 0.0_RP
    ATMOS_PHY_CP_RHOQV_t (:,:,:)   = 0.0_RP
    ATMOS_PHY_CP_RHOHYD_t(:,:,:,:) = 0.0_RP

    allocate( ATMOS_PHY_CP_w0mean        (KA,IA,JA) )
    allocate( ATMOS_PHY_CP_kf_nca        (IA,JA)    )
    ATMOS_PHY_CP_w0mean        (:,:,:) =    0.0_RP
    ATMOS_PHY_CP_kf_nca        (:,:)   = -100.0_RP

    ! for tendency restart
    VMAX_t = 3 + N_HYD
    allocate( VAR_t_NAME(VMAX_t) )
    allocate( VAR_t_DESC(VMAX_t) )
    allocate( VAR_t_UNIT(VMAX_t) )
    allocate( VAR_t_ID  (VMAX_t) )

    VAR_t_NAME(I_cp_dens_t) = 'DENS_t_CP'
    VAR_t_DESC(I_cp_dens_t) = 'tendency DENS in CP'
    VAR_t_UNIT(I_cp_dens_t) = 'kg/m3/s'
    VAR_t_NAME(I_cp_rhot_t) = 'RHOT_t_CP'
    VAR_t_DESC(I_cp_rhot_t) = 'tendency RHOT in CP'
    VAR_t_UNIT(I_cp_rhot_t) = 'K*kg/m3/s'

    VAR_t_NAME(I_cp_qv_t) = 'QV_t_CP'
    VAR_t_DESC(I_cp_qv_t) = 'tendency rho*QV in CP'
    VAR_t_UNIT(I_cp_qv_t) = 'kg/m3/s'
    do iq = 1, N_HYD
       VAR_t_NAME(3+iq) = trim(HYD_NAME(iq))//'_t_CP'
       VAR_t_DESC(3+iq) = 'tendency rho*'//trim(HYD_NAME(iq))//' in CP'
       VAR_t_UNIT(3+iq) = 'kg/m3/s'
    enddo


    allocate( ATMOS_PHY_CP_MFLX_cloudbase(IA,JA)    )
    allocate( ATMOS_PHY_CP_SFLX_rain     (IA,JA)    )
    allocate( ATMOS_PHY_CP_SFLX_snow     (IA,JA)    )
    allocate( ATMOS_PHY_CP_SFLX_ENGI     (IA,JA)    )
    allocate( ATMOS_PHY_CP_cloudtop      (IA,JA)    )
    allocate( ATMOS_PHY_CP_cloudbase     (IA,JA)    )
    allocate( ATMOS_PHY_CP_cldfrac_dp    (KA,IA,JA) )
    allocate( ATMOS_PHY_CP_cldfrac_sh    (KA,IA,JA) )
    ATMOS_PHY_CP_MFLX_cloudbase(:,:)   =    0.0_RP
    ATMOS_PHY_CP_SFLX_rain     (:,:)   =    0.0_RP
    ATMOS_PHY_CP_SFLX_snow     (:,:)   =    0.0_RP
    ATMOS_PHY_CP_SFLX_ENGI     (:,:)   =    0.0_RP
    ATMOS_PHY_CP_cloudtop      (:,:)   =    0.0_RP
    ATMOS_PHY_CP_cloudbase     (:,:)   =    0.0_RP
    ATMOS_PHY_CP_cldfrac_dp    (:,:,:) =    0.0_RP
    ATMOS_PHY_CP_cldfrac_sh    (:,:,:) =    0.0_RP


    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_CP_VARS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_PHY_CP_vars_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_CP_vars_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_CP_VARS. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_CP_VARS)

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_CP_vars_setup",*) '[ATMOS_PHY_CP] prognostic/diagnostic variables'
    LOG_INFO_CONT('(1x,A,A24,A,A48,A,A12,A)') &
               '      |', 'VARNAME                 ','|', &
               'DESCRIPTION                                     ', '[', 'UNIT        ', ']'
    do iv = 1, VMAX
       LOG_INFO_CONT('(1x,A,I3,A,A24,A,A48,A,A12,A)') &
                  'NO.',iv,'|',VAR_NAME(iv),'|',VAR_DESC(iv),'[',VAR_UNIT(iv),']'
    enddo

    ! tendency
    do iv = 1, VMAX_t
       LOG_INFO_CONT('(1x,A,I3,A,A24,A,A48,A,A12,A)') &
                  'NO.',iv+VMAX,'|',VAR_t_NAME(iv),'|',VAR_t_DESC(iv),'[',VAR_t_UNIT(iv),']'
    enddo

    LOG_NEWLINE
    if ( ATMOS_PHY_CP_RESTART_IN_BASENAME /= '' ) then
       LOG_INFO("ATMOS_PHY_CP_vars_setup",*) 'Restart input?  : YES, file = ', trim(ATMOS_PHY_CP_RESTART_IN_BASENAME)
       LOG_INFO("ATMOS_PHY_CP_vars_setup",*) 'Add timelabel?  : ', ATMOS_PHY_CP_RESTART_IN_POSTFIX_TIMELABEL
    else
       LOG_INFO("ATMOS_PHY_CP_vars_setup",*) 'Restart input?  : NO'
    endif
    if (       ATMOS_PHY_CP_RESTART_OUTPUT             &
         .AND. ATMOS_PHY_CP_RESTART_OUT_BASENAME /= '' ) then
       LOG_INFO("ATMOS_PHY_CP_vars_setup",*) 'Restart output? : YES, file = ', trim(ATMOS_PHY_CP_RESTART_OUT_BASENAME)
       LOG_INFO("ATMOS_PHY_CP_vars_setup",*) 'Add timelabel?  : ', ATMOS_PHY_CP_RESTART_OUT_POSTFIX_TIMELABEL
    else
       LOG_INFO("ATMOS_PHY_CP_vars_setup",*) 'Restart output? : NO'
       ATMOS_PHY_CP_RESTART_OUTPUT = .false.
    endif

    return
  end subroutine ATMOS_PHY_CP_vars_setup

  !-----------------------------------------------------------------------------
  !> HALO Communication
  subroutine ATMOS_PHY_CP_vars_fillhalo
    use scale_comm_cartesC, only: &
       COMM_vars8, &
       COMM_wait
    use scale_atmos_hydrometeor, only: &
       N_HYD
    implicit none

    integer :: i, j
    integer :: iq
    !---------------------------------------------------------------------------

    !$omp parallel do
    do j = JS, JE
    do i = IS, IE
       ATMOS_PHY_CP_w0mean (   1:KS-1,i,j) = ATMOS_PHY_CP_w0mean (KS,i,j)
       ATMOS_PHY_CP_w0mean (KE+1:KA  ,i,j) = ATMOS_PHY_CP_w0mean (KE,i,j)
       ATMOS_PHY_CP_DENS_t (   1:KS-1,i,j) = ATMOS_PHY_CP_DENS_t (KS,i,j)
       ATMOS_PHY_CP_DENS_t (KE+1:KA  ,i,j) = ATMOS_PHY_CP_DENS_t (KE,i,j)
       ATMOS_PHY_CP_RHOT_t (   1:KS-1,i,j) = ATMOS_PHY_CP_RHOT_t (KS,i,j)
       ATMOS_PHY_CP_RHOT_t (KE+1:KA  ,i,j) = ATMOS_PHY_CP_RHOT_t (KE,i,j)
       ATMOS_PHY_CP_RHOQV_t(   1:KS-1,i,j) = ATMOS_PHY_CP_RHOQV_t(KS,i,j)
       ATMOS_PHY_CP_RHOQV_t(KE+1:KA  ,i,j) = ATMOS_PHY_CP_RHOQV_t(KE,i,j)
    end do
    end do
    do iq = 1, N_HYD
       !$omp parallel do
       do j = JS, JE
       do i = IS, IE
          ATMOS_PHY_CP_RHOHYD_t(   1:KS-1,i,j,iq) = ATMOS_PHY_CP_RHOHYD_t(KS,i,j,iq)
          ATMOS_PHY_CP_RHOHYD_t(KE+1:KA  ,i,j,iq) = ATMOS_PHY_CP_RHOHYD_t(KE,i,j,iq)
       enddo
       enddo
    end do

    call COMM_vars8( ATMOS_PHY_CP_w0mean(:,:,:), 1 )
    call COMM_vars8( ATMOS_PHY_CP_kf_nca(:,:),   2 )

    ! tendency
    call COMM_vars8( ATMOS_PHY_CP_DENS_t (:,:,:), VMAX+1 )
    call COMM_vars8( ATMOS_PHY_CP_RHOT_t (:,:,:), VMAX+2 )
    call COMM_vars8( ATMOS_PHY_CP_RHOQV_t(:,:,:), VMAX+3 )
    do iq = 1, N_HYD
       call COMM_vars8( ATMOS_PHY_CP_RHOHYD_t(:,:,:,iq), VMAX+3+iq )
    enddo

    call COMM_wait ( ATMOS_PHY_CP_w0mean(:,:,:), 1 )
    call COMM_wait ( ATMOS_PHY_CP_kf_nca(:,:),   2 )

    call COMM_wait ( ATMOS_PHY_CP_DENS_t (:,:,:), VMAX+1 )
    call COMM_wait ( ATMOS_PHY_CP_RHOT_t (:,:,:), VMAX+2 )
    call COMM_wait ( ATMOS_PHY_CP_RHOQV_t(:,:,:), VMAX+3 )
    do iq = 1, N_HYD
       call COMM_wait ( ATMOS_PHY_CP_RHOHYD_t(:,:,:,iq), VMAX+3+iq )
    enddo

    return
  end subroutine ATMOS_PHY_CP_vars_fillhalo

  !-----------------------------------------------------------------------------
  !> Open restart file for read
  subroutine ATMOS_PHY_CP_vars_restart_open
    use scale_time, only: &
       TIME_gettimelabel
    use scale_file_cartesC, only: &
       FILE_CARTESC_open
    implicit none

    character(len=19)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_CP_vars_restart_open",*) 'Open restart file (ATMOS_PHY_CP) '

    if ( ATMOS_PHY_CP_RESTART_IN_BASENAME /= '' ) then

       if ( ATMOS_PHY_CP_RESTART_IN_POSTFIX_TIMELABEL ) then
          call TIME_gettimelabel( timelabel )
          basename = trim(ATMOS_PHY_CP_RESTART_IN_BASENAME)//'_'//trim(timelabel)
       else
          basename = trim(ATMOS_PHY_CP_RESTART_IN_BASENAME)
       endif

       LOG_INFO("ATMOS_PHY_CP_vars_restart_open",*) 'basename: ', trim(basename)

       call FILE_CARTESC_open( basename, restart_fid, aggregate=ATMOS_PHY_CP_RESTART_IN_AGGREGATE )
    else
       LOG_INFO("ATMOS_PHY_CP_vars_restart_open",*) 'restart file for ATMOS_PHY_CP is not specified.'
    endif

    return
  end subroutine ATMOS_PHY_CP_vars_restart_open

  !-----------------------------------------------------------------------------
  !> Read restart
  subroutine ATMOS_PHY_CP_vars_restart_read
    use scale_file, only: &
       FILE_get_AGGREGATE
    use scale_file_cartesC, only: &
       FILE_CARTESC_read, &
       FILE_CARTESC_flush
    use scale_atmos_hydrometeor, only: &
       N_HYD
    implicit none

    integer  :: i, j, iq
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then
       LOG_NEWLINE
       LOG_INFO("ATMOS_PHY_CP_vars_restart_read",*) 'Read from restart file (ATMOS_PHY_CP) '

       call FILE_CARTESC_read( restart_fid, VAR_NAME(1), 'ZXY', & ! [IN]
                               ATMOS_PHY_CP_w0mean(:,:,:)       ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_NAME(2), 'XY',  & ! [IN]
                               ATMOS_PHY_CP_kf_nca(:,:)         ) ! [OUT]
       ! tendency
       call FILE_CARTESC_read( restart_fid, VAR_t_NAME(1), 'ZXY', & ! [IN]
                               ATMOS_PHY_CP_DENS_t(:,:,:)         ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_t_NAME(2), 'ZXY', & ! [IN]
                               ATMOS_PHY_CP_RHOT_t(:,:,:)         ) ! [OUT]
       call FILE_CARTESC_read( restart_fid, VAR_t_NAME(3), 'ZXY', & ! [IN]
                               ATMOS_PHY_CP_RHOQV_t(:,:,:)        ) ! [OUT]
       do iq = 1, N_HYD
          call FILE_CARTESC_read( restart_fid, VAR_t_NAME(3+iq), 'ZXY', & ! [IN]
                                  ATMOS_PHY_CP_RHOHYD_t(:,:,:,iq)       ) ! [OUT]
       enddo

       if ( FILE_get_AGGREGATE(restart_fid) ) then
          call FILE_CARTESC_flush( restart_fid ) ! X/Y halos have been read from file

          ! fill K halos
          !$omp parallel do
          do j  = 1, JA
          do i  = 1, IA
             ATMOS_PHY_CP_w0mean (   1:KS-1,i,j) = ATMOS_PHY_CP_w0mean (KS,i,j)
             ATMOS_PHY_CP_w0mean (KE+1:KA  ,i,j) = ATMOS_PHY_CP_w0mean (KE,i,j)
             ATMOS_PHY_CP_DENS_t (   1:KS-1,i,j) = ATMOS_PHY_CP_DENS_t (KS,i,j)
             ATMOS_PHY_CP_DENS_t (KE+1:KA  ,i,j) = ATMOS_PHY_CP_DENS_t (KE,i,j)
             ATMOS_PHY_CP_RHOT_t (   1:KS-1,i,j) = ATMOS_PHY_CP_RHOT_t (KS,i,j)
             ATMOS_PHY_CP_RHOT_t (KE+1:KA  ,i,j) = ATMOS_PHY_CP_RHOT_t (KE,i,j)
             ATMOS_PHY_CP_RHOQV_t(   1:KS-1,i,j) = ATMOS_PHY_CP_RHOQV_t(KS,i,j)
             ATMOS_PHY_CP_RHOQV_t(KE+1:KA  ,i,j) = ATMOS_PHY_CP_RHOQV_t(KE,i,j)
          end do
          end do
          do iq = 1, N_HYD
             !$omp parallel do
             do j  = 1, JA
             do i  = 1, IA
                ATMOS_PHY_CP_RHOHYD_t(   1:KS-1,i,j,iq) = ATMOS_PHY_CP_RHOHYD_t(KS,i,j,iq)
                ATMOS_PHY_CP_RHOHYD_t(KE+1:KA,  i,j,iq) = ATMOS_PHY_CP_RHOHYD_t(KE,i,j,iq)
             enddo
             enddo
          enddo
       else
          call ATMOS_PHY_CP_vars_fillhalo
       end if

       call ATMOS_PHY_CP_vars_check

    else
       LOG_INFO("ATMOS_PHY_CP_vars_restart_read",*) 'invalid restart file ID for ATMOS_PHY_CP.'
    endif

    return
  end subroutine ATMOS_PHY_CP_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Create restart file
  subroutine ATMOS_PHY_CP_vars_restart_create
    use scale_time, only: &
       TIME_gettimelabel
    use scale_file_cartesC, only: &
       FILE_CARTESC_create
    implicit none

    character(len=19)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    if ( ATMOS_PHY_CP_RESTART_OUT_BASENAME /= '' ) then

       LOG_NEWLINE
       LOG_INFO("ATMOS_PHY_CP_vars_restart_create",*) 'Create restart file (ATMOS_PHY_AE) '

       if ( ATMOS_PHY_CP_RESTART_OUT_POSTFIX_TIMELABEL ) then
          call TIME_gettimelabel( timelabel )
          basename = trim(ATMOS_PHY_CP_RESTART_OUT_BASENAME)//'_'//trim(timelabel)
       else
          basename = trim(ATMOS_PHY_CP_RESTART_OUT_BASENAME)
       endif

       LOG_INFO("ATMOS_PHY_CP_vars_restart_create",*) 'basename: ', trim(basename)

       call FILE_CARTESC_create( &
            basename, ATMOS_PHY_CP_RESTART_OUT_TITLE, ATMOS_PHY_CP_RESTART_OUT_DTYPE, & ! [IN]
            restart_fid,                                                              & ! [OUT]
            aggregate=ATMOS_PHY_CP_RESTART_OUT_AGGREGATE                              ) ! [IN]

    endif

    return
  end subroutine ATMOS_PHY_CP_vars_restart_create

  !-----------------------------------------------------------------------------
  !> Exit netCDF define mode
  subroutine ATMOS_PHY_CP_vars_restart_enddef
    use scale_file_cartesC, only: &
       FILE_CARTESC_enddef
    implicit none

    if ( restart_fid /= -1 ) then
       call FILE_CARTESC_enddef( restart_fid ) ! [IN]
    endif

    return
  end subroutine ATMOS_PHY_CP_vars_restart_enddef

  !-----------------------------------------------------------------------------
  !> Close restart file
  subroutine ATMOS_PHY_CP_vars_restart_close
    use scale_file_cartesC, only: &
       FILE_CARTESC_close
    implicit none
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then
       LOG_NEWLINE
       LOG_INFO("ATMOS_PHY_CP_vars_restart_close",*) 'Close restart file (ATMOS_PHY_CP) '

       call FILE_CARTESC_close( restart_fid ) ! [IN]

       restart_fid = -1
    endif

    return
  end subroutine ATMOS_PHY_CP_vars_restart_close

  !-----------------------------------------------------------------------------
  !> Write restart
  subroutine ATMOS_PHY_CP_vars_restart_def_var
    use scale_file_cartesC, only: &
       FILE_CARTESC_def_var
    use scale_atmos_hydrometeor, only: &
       N_HYD
    implicit none

    integer :: i, iq
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then

       do i = 1, VMAX
          call FILE_CARTESC_def_var( restart_fid,          & ! [IN]
               VAR_NAME(i), VAR_DESC(i), VAR_UNIT(i),      & ! [IN]
               VAR_DIM(i), ATMOS_PHY_CP_RESTART_OUT_DTYPE, & ! [IN]
               VAR_ID(i)                                   ) ! [OUT]
       end do

       do i = 1, 3
          call FILE_CARTESC_def_var( restart_fid,           & ! [IN]
               VAR_t_NAME(i), VAR_t_DESC(i), VAR_t_UNIT(i), & ! [IN]
               'ZXY',  ATMOS_PHY_CP_RESTART_OUT_DTYPE,       & ! [IN]
               VAR_t_ID(i)                                  ) ! [OUT]
       end do

       do iq = 1, N_HYD
          call FILE_CARTESC_def_var( restart_fid,                    & ! [IN]
               VAR_t_NAME(3+iq), VAR_t_DESC(3+iq), VAR_t_UNIT(3+iq), & ! [IN]
               'ZXY',  ATMOS_PHY_CP_RESTART_OUT_DTYPE,               & ! [IN]
               VAR_t_ID(3+iq)                                        ) ! [OUT]
       enddo

    endif

    return
  end subroutine ATMOS_PHY_CP_vars_restart_def_var

  !-----------------------------------------------------------------------------
  !> Write restart
  subroutine ATMOS_PHY_CP_vars_restart_write
    use scale_file_cartesC, only: &
       FILE_CARTESC_write => FILE_CARTESC_write_var
    use scale_atmos_hydrometeor, only: &
       N_HYD
    implicit none

    integer  :: iq
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then

       call ATMOS_PHY_CP_vars_fillhalo

       call ATMOS_PHY_CP_vars_check

       call FILE_CARTESC_write( restart_fid, VAR_ID(1), ATMOS_PHY_CP_w0mean(:,:,:), & ! [IN]
                                VAR_NAME(1), 'ZXY' ) ! [IN]
       call FILE_CARTESC_write( restart_fid, VAR_ID(2), ATMOS_PHY_CP_kf_nca(:,:),   & ! [IN]
                                VAR_NAME(2), 'XY'  ) ! [IN]

       ! tendency
       call FILE_CARTESC_write( restart_fid, VAR_t_ID(1), ATMOS_PHY_CP_DENS_t(:,:,:), & ! [IN]
                          VAR_t_NAME(1), 'ZXY' ) ! [IN]
       call FILE_CARTESC_write( restart_fid, VAR_t_ID(2), ATMOS_PHY_CP_RHOT_t(:,:,:), & ! [IN]
                          VAR_t_NAME(2), 'ZXY' ) ! [IN]
       call FILE_CARTESC_write( restart_fid, VAR_t_ID(3), ATMOS_PHY_CP_RHOQV_t(:,:,:), & ! [IN]
                          VAR_t_NAME(3), 'ZXY' ) ! [IN]
       do iq = 1, N_HYD
          call FILE_CARTESC_write( restart_fid, VAR_t_ID(3+iq), ATMOS_PHY_CP_RHOHYD_t(:,:,:,iq), & ! [IN]
                             VAR_t_NAME(3+iq), 'ZXY' ) ! [IN]
       enddo

    endif

    return
  end subroutine ATMOS_PHY_CP_vars_restart_write

  subroutine ATMOS_PHY_CP_vars_check
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STATISTICS_total
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_AREA,    &
       ATMOS_GRID_CARTESC_REAL_TOTAREA, &
       ATMOS_GRID_CARTESC_REAL_VOL,     &
       ATMOS_GRID_CARTESC_REAL_TOTVOL
    use scale_atmos_hydrometeor, only: &
       N_HYD
    implicit none
    integer :: iq
    !---------------------------------------------------------------------------

    call VALCHECK( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                   ATMOS_PHY_CP_w0mean(:,:,:),         & ! (in)
                   -100.0_RP, 100.0_RP, VAR_NAME(1),   & ! (in)
                   __FILE__, __LINE__                  ) ! (in)
    call VALCHECK( IA, IS, IE, JA, JS, JE, &
                   ATMOS_PHY_CP_kf_nca(:,:),           & ! (in)
                   -100.0_RP, 1.0E5_RP, VAR_NAME(2),   & ! (in)
                   __FILE__, __LINE__                  ) ! (in)
    call VALCHECK( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                   ATMOS_PHY_CP_DENS_t(:,:,:),         & ! (in)
                   -1.0E0_RP, 1.0E0_RP, VAR_t_NAME(1), & ! (in)
                   __FILE__, __LINE__                  ) ! (in)
    call VALCHECK( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                   ATMOS_PHY_CP_RHOT_t(:,:,:),         & ! (in)
                   -1.0E3_RP, 1.0E3_RP, VAR_t_NAME(2), & ! (in)
                   __FILE__, __LINE__                  ) ! (in)
    call VALCHECK( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                   ATMOS_PHY_CP_RHOQV_t(:,:,:),        & ! (in)
                   -1.0E0_RP, 1.0E0_RP, VAR_t_NAME(3), & ! (in)
                   __FILE__, __LINE__                  ) ! (in)
    do iq = 1, N_HYD
       call VALCHECK( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                      ATMOS_PHY_CP_RHOHYD_t(:,:,:,iq),       & ! (in)
                      -1.0E0_RP, 1.0E0_RP, VAR_t_NAME(3+iq), & ! (in)
                      __FILE__, __LINE__                     ) ! (in)
    end do

    call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                           ATMOS_PHY_CP_w0mean        (:,:,:), VAR_NAME(1), & ! (in)
                           ATMOS_GRID_CARTESC_REAL_VOL(:,:,:),              & ! (in)
                           ATMOS_GRID_CARTESC_REAL_TOTVOL                   ) ! (in)
    call STATISTICS_total( IA, IS, IE, JA, JS, JE, &
                           ATMOS_PHY_CP_kf_nca        (:,:)  , VAR_NAME(2), & ! (in)
                           ATMOS_GRID_CARTESC_REAL_AREA(:,:),               & ! (in)
                           ATMOS_GRID_CARTESC_REAL_TOTAREA                  ) ! (in)
    ! tendency
    call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                           ATMOS_PHY_CP_DENS_t        (:,:,:), VAR_t_NAME(1), & ! (in)
                           ATMOS_GRID_CARTESC_REAL_VOL(:,:,:),                & ! (in)
                           ATMOS_GRID_CARTESC_REAL_TOTVOL                     ) ! (in)
    call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                           ATMOS_PHY_CP_RHOT_t        (:,:,:), VAR_t_NAME(2), & ! (in)
                           ATMOS_GRID_CARTESC_REAL_VOL(:,:,:),                & ! (in)
                           ATMOS_GRID_CARTESC_REAL_TOTVOL                     ) ! (in)
    call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                           ATMOS_PHY_CP_RHOQV_t       (:,:,:), VAR_t_NAME(3), & ! (in)
                           ATMOS_GRID_CARTESC_REAL_VOL(:,:,:),                & ! (in)
                           ATMOS_GRID_CARTESC_REAL_TOTVOL                     ) ! (in)
    do iq = 1, N_HYD
       call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                              ATMOS_PHY_CP_RHOHYD_t(:,:,:,iq), VAR_t_NAME(3+iq), & ! (in)
                              ATMOS_GRID_CARTESC_REAL_VOL(:,:,:),                & ! (in)
                              ATMOS_GRID_CARTESC_REAL_TOTVOL                     ) ! (in)
    enddo

    return
  end subroutine ATMOS_PHY_CP_vars_check

end module mod_atmos_phy_cp_vars
