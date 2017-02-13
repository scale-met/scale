!-------------------------------------------------------------------------------
!> module Atmosphere / Physics Cumulus
!!
!! @par Description
!!          Container for mod_atmos_phy_cp
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-05-04 (H.Yashiro)    [new]
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module mod_atmos_phy_cp_vars
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
  logical,               public :: ATMOS_PHY_CP_RESTART_OUTPUT                = .false.                !< output restart file?

  character(len=H_LONG),  public :: ATMOS_PHY_CP_RESTART_IN_BASENAME           = ''                     !< Basename of the input  file
  logical,                public :: ATMOS_PHY_CP_RESTART_IN_POSTFIX_TIMELABEL  = .false.                !< Add timelabel to the basename of input  file?
  character(len=H_LONG),  public :: ATMOS_PHY_CP_RESTART_OUT_BASENAME          = ''                     !< Basename of the output file
  logical,                public :: ATMOS_PHY_CP_RESTART_OUT_POSTFIX_TIMELABEL = .true.                 !< Add timelabel to the basename of output file?
  character(len=H_MID),   public :: ATMOS_PHY_CP_RESTART_OUT_TITLE             = 'ATMOS_PHY_CP restart' !< title    of the output file
  character(len=H_SHORT), public :: ATMOS_PHY_CP_RESTART_OUT_DTYPE             = 'DEFAULT'              !< REAL4 or REAL8

  real(RP), public, allocatable :: ATMOS_PHY_CP_DENS_t(:,:,:)    ! tendency DENS [kg/m3/s]
  real(RP), public, allocatable :: ATMOS_PHY_CP_MOMZ_t(:,:,:)    ! tendency MOMZ [kg/m2/s2]
  real(RP), public, allocatable :: ATMOS_PHY_CP_MOMX_t(:,:,:)    ! tendency MOMX [kg/m2/s2]
  real(RP), public, allocatable :: ATMOS_PHY_CP_MOMY_t(:,:,:)    ! tendency MOMY [kg/m2/s2]
  real(RP), public, allocatable :: ATMOS_PHY_CP_RHOT_t(:,:,:)    ! tendency RHOT [K*kg/m3/s]
  real(RP), public, allocatable :: ATMOS_PHY_CP_RHOQ_t(:,:,:,:)  ! tendency rho*QTRC [kg/kg/s]

  real(RP), public, allocatable :: ATMOS_PHY_CP_MFLX_cloudbase(:,:)   ! cloud base mass flux [kg/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_CP_SFLX_rain     (:,:)   ! convective rain [kg/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_CP_cloudtop      (:,:)   ! cloud top  height [m]
  real(RP), public, allocatable :: ATMOS_PHY_CP_cloudbase     (:,:)   ! cloud base height [m]
  real(RP), public, allocatable :: ATMOS_PHY_CP_cldfrac_dp    (:,:,:) ! cloud fraction (deep    convection) [0-1]
  real(RP), public, allocatable :: ATMOS_PHY_CP_cldfrac_sh    (:,:,:) ! cloud fraction (shallow convection) [0-1]
  ! only for K-F scheme
  real(RP), public, allocatable :: ATMOS_PHY_CP_kf_nca        (:,:)   ! advection/cumulus convection timescale/dt for KF[step]
  real(RP), public, allocatable :: ATMOS_PHY_CP_kf_w0avg      (:,:,:) ! running mean vertical wind velocity for KF[m/s]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,                private, parameter :: VMAX = 8       !< number of the variables
  integer,                private, parameter :: I_MFLX_cloudbase = 1
  integer,                private, parameter :: I_SFLX_convrain  = 2
  integer,                private, parameter :: I_cloudtop       = 3
  integer,                private, parameter :: I_cloudbase      = 4
  integer,                private, parameter :: I_cldfrac_dp     = 5
  integer,                private, parameter :: I_cldfrac_sh     = 6
  integer,                private, parameter :: I_kf_nca         = 7
  integer,                private, parameter :: I_kf_w0avg       = 8

  character(len=H_SHORT), private            :: VAR_NAME(VMAX) !< name  of the variables
  character(len=H_MID),   private            :: VAR_DESC(VMAX) !< desc. of the variables
  character(len=H_SHORT), private            :: VAR_UNIT(VMAX) !< unit  of the variables
  integer,                private            :: VAR_ID  (VMAX) !< ID    of the variables
  integer,                private            :: restart_fid = -1  ! file ID

  data VAR_NAME / 'MFLX_cloudbase',  &
                  'SFLX_convrain',   &
                  'cloudtop',        &
                  'cloudbase',       &
                  'cldfrac_dp',      &
                  'cldfrac_sh',      &
                  'kf_nca',          &
                  'kf_w0avg'         /
  data VAR_DESC / 'cloud base mass flux',                             &
                  'convective rain',                                  &
                  'cloud top height',                                 &
                  'cloud base height',                                &
                  'cloud fraction (deep convection)',                 &
                  'cloud fraction (shallow convection)',              &
                  'advection/cumulus convection timescale/dt for KF', &
                  'running mean vertical wind velocity for KF'        /
  data VAR_UNIT / 'kg/m2/s', &
                  'kg/m2/s', &
                  'm',       &
                  'm',       &
                  '0-1',     &
                  '0-1',     &
                  'step',    &
                  'm/s'      /

  ! tendency names
  integer,                private              :: VMAX_t       !< number of the tendency variables dens+rhot+QA_MP
  integer,                private              :: I_cp_dens_t = 1
  integer,                private              :: I_cp_rhot_t = 2

  character(len=H_SHORT), private, allocatable :: VAR_t_NAME(:) !< name  of the variables
  character(len=H_MID),   private, allocatable :: VAR_t_DESC(:) !< desc. of the variables
  character(len=H_SHORT), private, allocatable :: VAR_t_UNIT(:) !< unit  of the variables
  integer,                private, allocatable :: VAR_t_ID  (:) !< ID    of the variables

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_CP_vars_setup
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_atmos_phy_mp, only: &
       AQ_NAME => ATMOS_PHY_MP_NAME, &
       QA_MP
    implicit none

    NAMELIST / PARAM_ATMOS_PHY_CP_VARS / &
       ATMOS_PHY_CP_RESTART_IN_BASENAME,           &
       ATMOS_PHY_CP_RESTART_IN_POSTFIX_TIMELABEL,  &
       ATMOS_PHY_CP_RESTART_OUTPUT,                &
       ATMOS_PHY_CP_RESTART_OUT_BASENAME,          &
       ATMOS_PHY_CP_RESTART_OUT_POSTFIX_TIMELABEL, &
       ATMOS_PHY_CP_RESTART_OUT_TITLE,             &
       ATMOS_PHY_CP_RESTART_OUT_DTYPE

    integer :: ierr
    integer :: iv
    integer :: iq
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[VARS] / Categ[ATMOS PHY_CP] / Origin[SCALE-RM]'

    allocate( ATMOS_PHY_CP_DENS_t(KA,IA,JA)       )
    allocate( ATMOS_PHY_CP_MOMZ_t(KA,IA,JA)       )
    allocate( ATMOS_PHY_CP_MOMX_t(KA,IA,JA)       )
    allocate( ATMOS_PHY_CP_MOMY_t(KA,IA,JA)       )
    allocate( ATMOS_PHY_CP_RHOT_t(KA,IA,JA)       )
    allocate( ATMOS_PHY_CP_RHOQ_t(KA,IA,JA,QS_MP:QE_MP) )
    ATMOS_PHY_CP_DENS_t(:,:,:)   = 0.0_RP
    ATMOS_PHY_CP_MOMZ_t(:,:,:)   = UNDEF
    ATMOS_PHY_CP_MOMX_t(:,:,:)   = UNDEF
    ATMOS_PHY_CP_MOMY_t(:,:,:)   = UNDEF
    ATMOS_PHY_CP_RHOT_t(:,:,:)   = 0.0_RP
    ATMOS_PHY_CP_RHOQ_t(:,:,:,:) = 0.0_RP

    allocate( ATMOS_PHY_CP_MFLX_cloudbase(IA,JA)    )
    allocate( ATMOS_PHY_CP_SFLX_rain     (IA,JA)    )
    allocate( ATMOS_PHY_CP_cloudtop      (IA,JA)    )
    allocate( ATMOS_PHY_CP_cloudbase     (IA,JA)    )
    allocate( ATMOS_PHY_CP_cldfrac_dp    (KA,IA,JA) )
    allocate( ATMOS_PHY_CP_cldfrac_sh    (KA,IA,JA) )
    allocate( ATMOS_PHY_CP_kf_nca        (IA,JA)    )
    allocate( ATMOS_PHY_CP_kf_w0avg      (KA,IA,JA) )
    ATMOS_PHY_CP_MFLX_cloudbase(:,:)   =    0.0_RP
    ATMOS_PHY_CP_SFLX_rain     (:,:)   =    0.0_RP
    ATMOS_PHY_CP_cloudtop      (:,:)   =    0.0_RP
    ATMOS_PHY_CP_cloudbase     (:,:)   =    0.0_RP
    ATMOS_PHY_CP_cldfrac_dp    (:,:,:) =    0.0_RP
    ATMOS_PHY_CP_cldfrac_sh    (:,:,:) =    0.0_RP
    ATMOS_PHY_CP_kf_nca        (:,:)   = -100.0_RP
    ATMOS_PHY_CP_kf_w0avg      (:,:,:) =    0.0_RP

    ! for tendency restart
    VMAX_t = 2 + QA_MP
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

    do iq = 1, QA_MP
       VAR_t_NAME(2+iq) = trim(AQ_NAME(iq))//'_t_CP'
       VAR_t_DESC(2+iq) = 'tendency rho*'//trim(AQ_NAME(iq))//' in CP'
       VAR_t_UNIT(2+iq) = 'kg/m3/s'
    enddo

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_CP_VARS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_CP_VARS. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_ATMOS_PHY_CP_VARS)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** [ATMOS_PHY_CP] prognostic/diagnostic variables'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A24,A,A48,A,A12,A)') &
               '***       |', 'VARNAME                 ','|', &
               'DESCRIPTION                                     ', '[', 'UNIT        ', ']'
    do iv = 1, VMAX
       if( IO_L ) write(IO_FID_LOG,'(1x,A,I3,A,A24,A,A48,A,A12,A)') &
                  '*** NO.',iv,'|',VAR_NAME(iv),'|',VAR_DESC(iv),'[',VAR_UNIT(iv),']'
    enddo

    ! tendency
    do iv = 1, VMAX_t
       if( IO_L ) write(IO_FID_LOG,'(1x,A,I3,A,A24,A,A48,A,A12,A)') &
                  '*** NO.',iv+VMAX,'|',VAR_t_NAME(iv),'|',VAR_t_DESC(iv),'[',VAR_t_UNIT(iv),']'
    enddo

    if( IO_L ) write(IO_FID_LOG,*)
    if ( ATMOS_PHY_CP_RESTART_IN_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart input?  : YES, file = ', trim(ATMOS_PHY_CP_RESTART_IN_BASENAME)
       if( IO_L ) write(IO_FID_LOG,*) '*** Add timelabel?  : ', ATMOS_PHY_CP_RESTART_IN_POSTFIX_TIMELABEL
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart input?  : NO'
    endif
    if (       ATMOS_PHY_CP_RESTART_OUTPUT             &
         .AND. ATMOS_PHY_CP_RESTART_OUT_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart output? : YES, file = ', trim(ATMOS_PHY_CP_RESTART_OUT_BASENAME)
       if( IO_L ) write(IO_FID_LOG,*) '*** Add timelabel?  : ', ATMOS_PHY_CP_RESTART_OUT_POSTFIX_TIMELABEL
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart output? : NO'
       ATMOS_PHY_CP_RESTART_OUTPUT = .false.
    endif

    return
  end subroutine ATMOS_PHY_CP_vars_setup

  !-----------------------------------------------------------------------------
  !> HALO Communication
  subroutine ATMOS_PHY_CP_vars_fillhalo
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    use scale_atmos_phy_mp, only: &
       QA_MP
    implicit none

    integer :: i, j
    integer :: iq
    !---------------------------------------------------------------------------

    do j  = JS, JE
    do i  = IS, IE
       ATMOS_PHY_CP_cldfrac_dp(   1:KS-1,i,j) = ATMOS_PHY_CP_cldfrac_dp(KS,i,j)
       ATMOS_PHY_CP_cldfrac_dp(KE+1:KA,  i,j) = ATMOS_PHY_CP_cldfrac_dp(KE,i,j)
       ATMOS_PHY_CP_cldfrac_sh(   1:KS-1,i,j) = ATMOS_PHY_CP_cldfrac_sh(KS,i,j)
       ATMOS_PHY_CP_cldfrac_sh(KE+1:KA,  i,j) = ATMOS_PHY_CP_cldfrac_sh(KE,i,j)
       ATMOS_PHY_CP_kf_w0avg  (   1:KS-1,i,j) = ATMOS_PHY_CP_kf_w0avg  (KS,i,j)
       ATMOS_PHY_CP_kf_w0avg  (KE+1:KA,  i,j) = ATMOS_PHY_CP_kf_w0avg  (KE,i,j)
       ATMOS_PHY_CP_DENS_t    (   1:KS-1,i,j) = ATMOS_PHY_CP_DENS_t    (KS,i,j)
       ATMOS_PHY_CP_DENS_t    (KE+1:KA  ,i,j) = ATMOS_PHY_CP_DENS_t    (KE,i,j)
       ATMOS_PHY_CP_RHOT_t    (   1:KS-1,i,j) = ATMOS_PHY_CP_RHOT_t    (KS,i,j)
       ATMOS_PHY_CP_RHOT_t    (KE+1:KA  ,i,j) = ATMOS_PHY_CP_RHOT_t    (KE,i,j)
    enddo
    enddo

    do iq = 1, QA_MP
       do j  = JS, JE
       do i  = IS, IE
          ATMOS_PHY_CP_RHOQ_t(   1:KS-1,i,j,iq) = ATMOS_PHY_CP_RHOQ_t(KS,i,j,iq)
          ATMOS_PHY_CP_RHOQ_t(KE+1:KA  ,i,j,iq) = ATMOS_PHY_CP_RHOQ_t(KE,i,j,iq)
       enddo
       enddo
    end do

    call COMM_vars8( ATMOS_PHY_CP_MFLX_cloudbase (:,:)  , 1 )
    call COMM_vars8( ATMOS_PHY_CP_SFLX_rain      (:,:)  , 2 )
    call COMM_vars8( ATMOS_PHY_CP_cloudtop       (:,:)  , 3 )
    call COMM_vars8( ATMOS_PHY_CP_cloudbase      (:,:)  , 4 )
    call COMM_vars8( ATMOS_PHY_CP_cldfrac_dp     (:,:,:), 5 )
    call COMM_vars8( ATMOS_PHY_CP_cldfrac_sh     (:,:,:), 6 )
    call COMM_vars8( ATMOS_PHY_CP_kf_nca         (:,:)  , 7 )
    call COMM_vars8( ATMOS_PHY_CP_kf_w0avg       (:,:,:), 8 )

    ! tendency
    call COMM_vars8( ATMOS_PHY_CP_DENS_t(:,:,:), VMAX+1 )
    call COMM_vars8( ATMOS_PHY_CP_RHOT_t(:,:,:), VMAX+2 )

    do iq = 1, QA_MP
       call COMM_vars8( ATMOS_PHY_CP_RHOQ_t(:,:,:,iq), VMAX+2+iq )
    enddo

    call COMM_wait ( ATMOS_PHY_CP_MFLX_cloudbase (:,:)  , 1 )
    call COMM_wait ( ATMOS_PHY_CP_SFLX_rain      (:,:)  , 2 )
    call COMM_wait ( ATMOS_PHY_CP_cloudtop       (:,:)  , 3 )
    call COMM_wait ( ATMOS_PHY_CP_cloudbase      (:,:)  , 4 )
    call COMM_wait ( ATMOS_PHY_CP_cldfrac_dp     (:,:,:), 5 )
    call COMM_wait ( ATMOS_PHY_CP_cldfrac_sh     (:,:,:), 6 )
    call COMM_wait ( ATMOS_PHY_CP_kf_nca         (:,:)  , 7 )
    call COMM_wait ( ATMOS_PHY_CP_kf_w0avg       (:,:,:), 8 )

    call COMM_wait ( ATMOS_PHY_CP_DENS_t(:,:,:), VMAX+1 )
    call COMM_wait ( ATMOS_PHY_CP_RHOT_t(:,:,:), VMAX+2 )

    do iq = 1, QA_MP
       call COMM_wait ( ATMOS_PHY_CP_RHOQ_t(:,:,:,iq), VMAX+2+iq )
    enddo

    return
  end subroutine ATMOS_PHY_CP_vars_fillhalo

  !-----------------------------------------------------------------------------
  !> Open restart file for read
  subroutine ATMOS_PHY_CP_vars_restart_open
    use scale_time, only: &
       TIME_gettimelabel
    use scale_fileio, only: &
       FILEIO_open
    implicit none

    character(len=19)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Open restart file (ATMOS_PHY_CP) ***'

    if ( ATMOS_PHY_CP_RESTART_IN_BASENAME /= '' ) then

       if ( ATMOS_PHY_CP_RESTART_IN_POSTFIX_TIMELABEL ) then
          call TIME_gettimelabel( timelabel )
          basename = trim(ATMOS_PHY_CP_RESTART_IN_BASENAME)//'_'//trim(timelabel)
       else
          basename = trim(ATMOS_PHY_CP_RESTART_IN_BASENAME)
       endif

       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(basename)

       call FILEIO_open( restart_fid, basename )
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** restart file for ATMOS_PHY_CP is not specified.'
    endif

    return
  end subroutine ATMOS_PHY_CP_vars_restart_open

  !-----------------------------------------------------------------------------
  !> Read restart
  subroutine ATMOS_PHY_CP_vars_restart_read
    use scale_rm_statistics, only: &
       STATISTICS_checktotal, &
       STAT_total
    use scale_fileio, only: &
       FILEIO_read, &
       FILEIO_flush
    use scale_atmos_phy_mp, only: &
       QA_MP
    implicit none

    real(RP) :: total
    integer  :: i, j, iq
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Read from restart file (ATMOS_PHY_CP) ***'

       call FILEIO_read( ATMOS_PHY_CP_MFLX_cloudbase(:,:),                            & ! [OUT]
                         restart_fid, VAR_NAME(1), 'XY',  step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_CP_SFLX_rain(:,:),                                 & ! [OUT]
                         restart_fid, VAR_NAME(2), 'XY',  step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_CP_cloudtop(:,:),                                  & ! [OUT]
                         restart_fid, VAR_NAME(3), 'XY',  step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_CP_cloudbase(:,:),                                 & ! [OUT]
                         restart_fid, VAR_NAME(4), 'XY',  step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_CP_cldfrac_dp(:,:,:),                              & ! [OUT]
                         restart_fid, VAR_NAME(5), 'ZXY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_CP_cldfrac_sh(:,:,:),                              & ! [OUT]
                         restart_fid, VAR_NAME(6), 'ZXY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_CP_kf_nca(:,:),                                    & ! [OUT]
                         restart_fid, VAR_NAME(7), 'XY',  step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_CP_kf_w0avg(:,:,:),                                & ! [OUT]
                         restart_fid, VAR_NAME(8), 'ZXY', step=1 ) ! [IN]
       ! tendency
       call FILEIO_read( ATMOS_PHY_CP_DENS_t(:,:,:),                                    & ! [OUT]
                         restart_fid, VAR_t_NAME(1), 'ZXY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_CP_RHOT_t(:,:,:),                                    & ! [OUT]
                         restart_fid, VAR_t_NAME(2), 'ZXY', step=1 ) ! [IN]
       do iq = 1, QA_MP
          call FILEIO_read( ATMOS_PHY_CP_RHOQ_t(:,:,:,iq),                                    & ! [OUT]
                            restart_fid, VAR_t_NAME(2+iq), 'ZXY', step=1 ) ! [IN]
       enddo

       if ( IO_AGGREGATE ) then
          call FILEIO_flush( restart_fid ) ! X/Y halos have been read from file

          ! fill K halos
          do j  = 1, JA
          do i  = 1, IA
             ATMOS_PHY_CP_cldfrac_dp(   1:KS-1,i,j) = ATMOS_PHY_CP_cldfrac_dp(KS,i,j)
             ATMOS_PHY_CP_cldfrac_sh(   1:KS-1,i,j) = ATMOS_PHY_CP_cldfrac_sh(KS,i,j)
             ATMOS_PHY_CP_kf_w0avg  (   1:KS-1,i,j) = ATMOS_PHY_CP_kf_w0avg  (KS,i,j)
             ATMOS_PHY_CP_DENS_t    (   1:KS-1,i,j) = ATMOS_PHY_CP_DENS_t    (KS,i,j)
             ATMOS_PHY_CP_RHOT_t    (   1:KS-1,i,j) = ATMOS_PHY_CP_RHOT_t    (KS,i,j)
             ATMOS_PHY_CP_cldfrac_dp(KE+1:KA,  i,j) = ATMOS_PHY_CP_cldfrac_dp(KE,i,j)
             ATMOS_PHY_CP_cldfrac_sh(KE+1:KA,  i,j) = ATMOS_PHY_CP_cldfrac_sh(KE,i,j)
             ATMOS_PHY_CP_kf_w0avg  (KE+1:KA,  i,j) = ATMOS_PHY_CP_kf_w0avg  (KE,i,j)
             ATMOS_PHY_CP_DENS_t    (KE+1:KA,  i,j) = ATMOS_PHY_CP_DENS_t    (KE,i,j)
             ATMOS_PHY_CP_RHOT_t    (KE+1:KA,  i,j) = ATMOS_PHY_CP_RHOT_t    (KE,i,j)
          enddo
          enddo

          do iq = 1, QA_MP
             do j  = 1, JA
             do i  = 1, IA
                ATMOS_PHY_CP_RHOQ_t(   1:KS-1,i,j,iq) = ATMOS_PHY_CP_RHOQ_t(KS,i,j,iq)
                ATMOS_PHY_CP_RHOQ_t(KE+1:KA,  i,j,iq) = ATMOS_PHY_CP_RHOQ_t(KE,i,j,iq)
             enddo
             enddo
          enddo
       else
          call ATMOS_PHY_CP_vars_fillhalo
       end if

       if ( STATISTICS_checktotal ) then
          call STAT_total( total, ATMOS_PHY_CP_MFLX_cloudbase(:,:)  , VAR_NAME(1) )
          call STAT_total( total, ATMOS_PHY_CP_SFLX_rain     (:,:)  , VAR_NAME(2) )
          call STAT_total( total, ATMOS_PHY_CP_cloudtop      (:,:)  , VAR_NAME(3) )
          call STAT_total( total, ATMOS_PHY_CP_cloudbase     (:,:)  , VAR_NAME(4) )
          call STAT_total( total, ATMOS_PHY_CP_cldfrac_dp    (:,:,:), VAR_NAME(5) )
          call STAT_total( total, ATMOS_PHY_CP_cldfrac_sh    (:,:,:), VAR_NAME(6) )
          call STAT_total( total, ATMOS_PHY_CP_kf_nca        (:,:)  , VAR_NAME(7) )
          call STAT_total( total, ATMOS_PHY_CP_kf_w0avg      (:,:,:), VAR_NAME(8) )
          ! tendency
          call STAT_total( total, ATMOS_PHY_CP_DENS_t        (:,:,:), VAR_t_NAME(1) )
          call STAT_total( total, ATMOS_PHY_CP_RHOT_t        (:,:,:), VAR_t_NAME(2) )
          do iq = 1, QA_MP
             call STAT_total( total, ATMOS_PHY_CP_RHOQ_t(:,:,:,iq), VAR_t_NAME(2+iq) )
          enddo
       endif
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** invalid restart file ID for ATMOS_PHY_CP.'
    endif

    return
  end subroutine ATMOS_PHY_CP_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Create restart file
  subroutine ATMOS_PHY_CP_vars_restart_create
    use scale_time, only: &
       TIME_gettimelabel
    use scale_fileio, only: &
       FILEIO_create
    implicit none

    character(len=19)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    if ( ATMOS_PHY_CP_RESTART_OUT_BASENAME /= '' ) then

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Create restart file (ATMOS_PHY_AE) ***'

       if ( ATMOS_PHY_CP_RESTART_OUT_POSTFIX_TIMELABEL ) then
          call TIME_gettimelabel( timelabel )
          basename = trim(ATMOS_PHY_CP_RESTART_OUT_BASENAME)//'_'//trim(timelabel)
       else
          basename = trim(ATMOS_PHY_CP_RESTART_OUT_BASENAME)
       endif

       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(basename)

       call FILEIO_create( restart_fid,                                                             & ! [OUT]
                           basename, ATMOS_PHY_CP_RESTART_OUT_TITLE, ATMOS_PHY_CP_RESTART_OUT_DTYPE ) ! [IN]

    endif

    return
  end subroutine ATMOS_PHY_CP_vars_restart_create

  !-----------------------------------------------------------------------------
  !> Exit netCDF define mode
  subroutine ATMOS_PHY_CP_vars_restart_enddef
    use scale_fileio, only: &
       FILEIO_enddef
    implicit none

    if ( restart_fid /= -1 ) then
       call FILEIO_enddef( restart_fid ) ! [IN]
    endif

    return
  end subroutine ATMOS_PHY_CP_vars_restart_enddef

  !-----------------------------------------------------------------------------
  !> Close restart file
  subroutine ATMOS_PHY_CP_vars_restart_close
    use scale_fileio, only: &
       FILEIO_close
    implicit none
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Close restart file (ATMOS_PHY_CP) ***'

       call FILEIO_close( restart_fid ) ! [IN]

       restart_fid = -1
    endif

    return
  end subroutine ATMOS_PHY_CP_vars_restart_close

  !-----------------------------------------------------------------------------
  !> Write restart
  subroutine ATMOS_PHY_CP_vars_restart_def_var
    use scale_fileio, only: &
       FILEIO_def_var
    use scale_atmos_phy_mp, only: &
       QA_MP
    implicit none

    integer :: iq
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then

       call FILEIO_def_var( restart_fid, VAR_ID(1), VAR_NAME(1), VAR_DESC(1),   &
                            VAR_UNIT(1), 'XY',  ATMOS_PHY_CP_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_def_var( restart_fid, VAR_ID(2), VAR_NAME(2), VAR_DESC(2),   &
                            VAR_UNIT(2), 'XY',  ATMOS_PHY_CP_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_def_var( restart_fid, VAR_ID(3), VAR_NAME(3), VAR_DESC(3),   &
                            VAR_UNIT(3), 'XY',  ATMOS_PHY_CP_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_def_var( restart_fid, VAR_ID(4), VAR_NAME(4), VAR_DESC(4),   &
                            VAR_UNIT(4), 'XY',  ATMOS_PHY_CP_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_def_var( restart_fid, VAR_ID(5), VAR_NAME(5), VAR_DESC(5),   &
                            VAR_UNIT(5), 'ZXY', ATMOS_PHY_CP_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_def_var( restart_fid, VAR_ID(6), VAR_NAME(6), VAR_DESC(6),   &
                            VAR_UNIT(6), 'ZXY', ATMOS_PHY_CP_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_def_var( restart_fid, VAR_ID(7), VAR_NAME(7), VAR_DESC(7),   &
                            VAR_UNIT(7), 'XY',  ATMOS_PHY_CP_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_def_var( restart_fid, VAR_ID(8), VAR_NAME(8), VAR_DESC(8),   &
                            VAR_UNIT(8), 'ZXY', ATMOS_PHY_CP_RESTART_OUT_DTYPE  ) ! [IN]

       call FILEIO_def_var( restart_fid, VAR_t_ID(1), VAR_t_NAME(1), VAR_t_DESC(1), &
                            VAR_t_UNIT(1), 'ZXY',  ATMOS_PHY_CP_RESTART_OUT_DTYPE   ) ! [IN]
       call FILEIO_def_var( restart_fid, VAR_t_ID(2), VAR_t_NAME(2), VAR_t_DESC(2), &
                            VAR_t_UNIT(2), 'ZXY',  ATMOS_PHY_CP_RESTART_OUT_DTYPE   ) ! [IN]

       do iq = 1, QA_MP
          call FILEIO_def_var( restart_fid, VAR_t_ID(2+iq), VAR_t_NAME(2+iq), VAR_t_DESC(2+iq), &
                               VAR_t_UNIT(2+iq), 'ZXY',  ATMOS_PHY_CP_RESTART_OUT_DTYPE         ) ! [IN]
       enddo

    endif

    return
  end subroutine ATMOS_PHY_CP_vars_restart_def_var

  !-----------------------------------------------------------------------------
  !> Write restart
  subroutine ATMOS_PHY_CP_vars_restart_write
    use scale_rm_statistics, only: &
       STATISTICS_checktotal, &
       STAT_total
    use scale_fileio, only: &
       FILEIO_write => FILEIO_write_var
    use scale_atmos_phy_mp, only: &
       QA_MP
    implicit none

    real(RP) :: total
    integer  :: iq
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then

       call ATMOS_PHY_CP_vars_fillhalo

       if ( STATISTICS_checktotal ) then
          call STAT_total( total, ATMOS_PHY_CP_MFLX_cloudbase(:,:)  , VAR_NAME(1) )
          call STAT_total( total, ATMOS_PHY_CP_SFLX_rain     (:,:)  , VAR_NAME(2) )
          call STAT_total( total, ATMOS_PHY_CP_cloudtop      (:,:)  , VAR_NAME(3) )
          call STAT_total( total, ATMOS_PHY_CP_cloudbase     (:,:)  , VAR_NAME(4) )
          call STAT_total( total, ATMOS_PHY_CP_cldfrac_dp    (:,:,:), VAR_NAME(5) )
          call STAT_total( total, ATMOS_PHY_CP_cldfrac_sh    (:,:,:), VAR_NAME(6) )
          call STAT_total( total, ATMOS_PHY_CP_kf_nca        (:,:)  , VAR_NAME(7) )
          call STAT_total( total, ATMOS_PHY_CP_kf_w0avg      (:,:,:), VAR_NAME(8) )
          ! tendency
          call STAT_total( total, ATMOS_PHY_CP_DENS_t        (:,:,:), VAR_t_NAME(1) )
          call STAT_total( total, ATMOS_PHY_CP_RHOT_t        (:,:,:), VAR_t_NAME(2) )
          do iq = 1, QA_MP
             call STAT_total( total, ATMOS_PHY_CP_RHOQ_t(:,:,:,iq), VAR_t_NAME(2+iq) )
          enddo
       endif

       call FILEIO_write( restart_fid, VAR_ID(1), ATMOS_PHY_CP_MFLX_cloudbase(:,:), & ! [IN]
                          VAR_NAME(1), 'XY' ) ! [IN]
       call FILEIO_write( restart_fid, VAR_ID(2), ATMOS_PHY_CP_SFLX_rain(:,:), & ! [IN]
                          VAR_NAME(2), 'XY' ) ! [IN]
       call FILEIO_write( restart_fid, VAR_ID(3), ATMOS_PHY_CP_cloudtop(:,:), & ! [IN]
                          VAR_NAME(3), 'XY' ) ! [IN]
       call FILEIO_write( restart_fid, VAR_ID(4), ATMOS_PHY_CP_cloudbase(:,:), & ! [IN]
                          VAR_NAME(4), 'XY' ) ! [IN]
       call FILEIO_write( restart_fid, VAR_ID(5), ATMOS_PHY_CP_cldfrac_dp(:,:,:), & ! [IN]
                          VAR_NAME(5), 'ZXY' ) ! [IN]
       call FILEIO_write( restart_fid, VAR_ID(6), ATMOS_PHY_CP_cldfrac_sh(:,:,:), & ! [IN]
                          VAR_NAME(6), 'ZXY' ) ! [IN]
       call FILEIO_write( restart_fid, VAR_ID(7), ATMOS_PHY_CP_kf_nca(:,:), & ! [IN]
                          VAR_NAME(7), 'XY' ) ! [IN]
       call FILEIO_write( restart_fid, VAR_ID(8), ATMOS_PHY_CP_kf_w0avg(:,:,:), & ! [IN]
                          VAR_NAME(8), 'ZXY' ) ! [IN]

       ! tendency
       call FILEIO_write( restart_fid, VAR_t_ID(1), ATMOS_PHY_CP_DENS_t(:,:,:), & ! [IN]
                          VAR_t_NAME(1), 'ZXY' ) ! [IN]
       call FILEIO_write( restart_fid, VAR_t_ID(2), ATMOS_PHY_CP_RHOT_t(:,:,:), & ! [IN]
                          VAR_t_NAME(2), 'ZXY' ) ! [IN]
       do iq = 1, QA_MP
          call FILEIO_write( restart_fid, VAR_t_ID(2+iq), ATMOS_PHY_CP_RHOQ_t(:,:,:,iq), & ! [IN]
                             VAR_t_NAME(2+iq), 'ZXY' ) ! [IN]
       enddo

    endif

    return
  end subroutine ATMOS_PHY_CP_vars_restart_write

end module mod_atmos_phy_cp_vars
