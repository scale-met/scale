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
  public :: ATMOS_PHY_CP_vars_restart_def_var
  public :: ATMOS_PHY_CP_vars_restart_enddef
  public :: ATMOS_PHY_CP_vars_restart_write_var
  public :: ATMOS_PHY_CP_vars_restart_close

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical,               public :: ATMOS_PHY_CP_RESTART_OUTPUT       = .false.                !< output restart file?

  character(len=H_LONG), public :: ATMOS_PHY_CP_RESTART_IN_BASENAME  = ''                     !< basename of the restart file
  character(len=H_LONG), public :: ATMOS_PHY_CP_RESTART_OUT_BASENAME = ''                     !< basename of the output file
  character(len=H_MID),  public :: ATMOS_PHY_CP_RESTART_OUT_TITLE    = 'ATMOS_PHY_CP restart' !< title    of the output file
  character(len=H_MID),  public :: ATMOS_PHY_CP_RESTART_OUT_DTYPE    = 'DEFAULT'              !< REAL4 or REAL8

  real(RP), public, allocatable :: ATMOS_PHY_CP_DENS_t(:,:,:)   ! tendency DENS [kg/m3/s]
  real(RP), public, allocatable :: ATMOS_PHY_CP_MOMZ_t(:,:,:)   ! tendency MOMZ [kg/m2/s2]
  real(RP), public, allocatable :: ATMOS_PHY_CP_MOMX_t(:,:,:)   ! tendency MOMX [kg/m2/s2]
  real(RP), public, allocatable :: ATMOS_PHY_CP_MOMY_t(:,:,:)   ! tendency MOMY [kg/m2/s2]
  real(RP), public, allocatable :: ATMOS_PHY_CP_RHOT_t(:,:,:)   ! tendency RHOT [K*kg/m3/s]
  real(RP), public, allocatable :: ATMOS_PHY_CP_RHOQ_t(:,:,:,:) ! tendency rho*QTRC [kg/kg/s]

  real(RP), public, allocatable :: ATMOS_PHY_CP_MFLX_cloudbase(:,:) ! cloud base mass flux [kg/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_CP_SFLX_convrain(:,:)  ! convective rain [kg/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_CP_cloudtop(:,:)       ! cloud top height [m]
  real(RP), public, allocatable :: ATMOS_PHY_CP_cloudbase(:,:)      ! cloud base height [m]
  real(RP), public, allocatable :: ATMOS_PHY_CP_cldfrac_dp(:,:,:)   ! cloud fraction (deep convection) [0-1]
  real(RP), public, allocatable :: ATMOS_PHY_CP_cldfrac_sh(:,:,:)   ! cloud fraction (shallow convection) [0-1]
  ! only use kf scheme
  real(RP), public, allocatable :: ATMOS_PHY_CP_kf_nca(:,:)         ! advection/cumulus convection timescale/dt for KF[step]
  real(RP), public, allocatable :: ATMOS_PHY_CP_kf_w0avg(:,:,:)       ! rannning mean vertical wind velocity for KF[m/s]



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

  data VAR_NAME / 'MFLX_cloudbase' , &
                  'SFLX_convrain'  , &
                  'cloudtop'       , &
                  'cloudbase'      , &
                  'cldfrac_dp'     , &
                  'cldfrac_sh'     , &
                  'kf_nca'         , &
                  'kf_w0avg'        &
                  /
  data VAR_DESC / 'cloud base mass flux', &
       'convective rain', &
       'cloud top height', &
       'cloud base height', &
       'cloud fraction (deep convection)', &
       'cloud fraction (shallow convection)', &
       'advection/cumulus convection timescale/dt for KF', &
       'rannning mean vertical wind velocity for KF' &
       /
  data VAR_UNIT / 'kg/m2/s', &
       'kg/m2/s', &
       'm', &
       'm', &
       '0-1', &
       '0-1', &
       'step', &
       'm/s' &
       /

  ! tendency names
  integer,                private, save      :: VMAX_tend = 2!+QA   !< number of the tendency variables dens+rhot+QA
  integer,                private, save      :: I_cp_dens_t  = 1
  integer,                private, save      :: I_cp_rhot_t  = 2
  integer,private,    allocatable, save      :: I_cp_rhoq_t(:)!(QA) = 2+I_qv,2+i_qc....
  character(len=H_SHORT), private,allocatable, save      :: VAR_tend_NAME(:)!(VMAX_tend) !< name  of the variables
  character(len=H_MID),   private,allocatable, save      :: VAR_tend_DESC(:)!(VMAX_tend) !< desc. of the variables
  character(len=H_SHORT), private,allocatable, save      :: VAR_tend_UNIT(:)!(VMAX_tend) !< unit  of the variables

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_CP_vars_setup
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    implicit none

    NAMELIST / PARAM_ATMOS_PHY_CP_VARS / &
       ATMOS_PHY_CP_RESTART_IN_BASENAME,  &
       ATMOS_PHY_CP_RESTART_OUTPUT,       &
       ATMOS_PHY_CP_RESTART_OUT_BASENAME, &
       ATMOS_PHY_CP_RESTART_OUT_TITLE,    &
       ATMOS_PHY_CP_RESTART_OUT_DTYPE

    integer :: ierr
    integer :: iv
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[VARS] / Categ[ATMOS PHY_CP] / Origin[SCALE-RM]'

    allocate( ATMOS_PHY_CP_DENS_t(KA,IA,JA)    )
    allocate( ATMOS_PHY_CP_MOMZ_t(KA,IA,JA)    )
    allocate( ATMOS_PHY_CP_MOMX_t(KA,IA,JA)    )
    allocate( ATMOS_PHY_CP_MOMY_t(KA,IA,JA)    )
    allocate( ATMOS_PHY_CP_RHOT_t(KA,IA,JA)    )
    allocate( ATMOS_PHY_CP_RHOQ_t(KA,IA,JA,QA) )
    ATMOS_PHY_CP_DENS_t(:,:,:)   = 0._RP!UNDEF
    ATMOS_PHY_CP_MOMZ_t(:,:,:)   = UNDEF
    ATMOS_PHY_CP_MOMX_t(:,:,:)   = UNDEF
    ATMOS_PHY_CP_MOMY_t(:,:,:)   = UNDEF
    ATMOS_PHY_CP_RHOT_t(:,:,:)   = 0._RP!UNDEF
    ATMOS_PHY_CP_RHOQ_t(:,:,:,:) = 0._RP!UNDEF

    allocate( ATMOS_PHY_CP_MFLX_cloudbase(IA,JA) )
    allocate( ATMOS_PHY_CP_SFLX_convrain (IA,JA) )
    allocate( ATMOS_PHY_CP_cloudtop      (IA,JA) )
    allocate( ATMOS_PHY_CP_cloudbase     (IA,JA) )
    allocate( ATMOS_PHY_CP_cldfrac_dp (KA,IA,JA) )
    allocate( ATMOS_PHY_CP_cldfrac_sh (KA,IA,JA) )
    allocate( ATMOS_PHY_CP_kf_nca        (IA,JA) )
    allocate( ATMOS_PHY_CP_kf_w0avg   (KA,IA,JA) )
    ATMOS_PHY_CP_MFLX_cloudbase (:,:)   = 0._RP
    ATMOS_PHY_CP_SFLX_convrain  (:,:)   = 0._RP
    ATMOS_PHY_CP_cloudtop       (:,:)   = 0._RP
    ATMOS_PHY_CP_cloudbase      (:,:)   = 0._RP
    ATMOS_PHY_CP_cldfrac_dp     (:,:,:) = 0._RP
    ATMOS_PHY_CP_cldfrac_sh     (:,:,:) = 0._RP
    ATMOS_PHY_CP_kf_nca         (:,:)   = -100._RP
    ATMOS_PHY_CP_kf_w0avg       (:,:,:) = 0._RP
    ! for tendency restart
    allocate( I_cp_rhoq_t(QA))
    I_cp_rhoq_t(:) = (/(I_cp_rhot_t+iv,iv=1,QA)/)
    VMAX_tend = VMAX_tend + QA
    allocate (VAR_tend_NAME(VMAX_tend))
    allocate (VAR_tend_DESC(VMAX_tend))
    allocate (VAR_tend_UNIT(VMAX_tend))
    VAR_tend_NAME(I_cp_dens_t) = 'DENS_t_CP' ; VAR_tend_DESC(I_cp_dens_t) = 'tendency DENS in CP'
    VAR_tend_UNIT(I_cp_dens_t) = 'kg/m3/s'
    VAR_tend_NAME(I_cp_rhot_t) = 'RHOT_t_CP' ; VAR_tend_DESC(I_cp_rhot_t) = 'tendency RHOT in CP'
    VAR_tend_UNIT(I_cp_rhot_t) = 'K*kg/m3/s'
    ! moisture
    do iv = 1,QA
       VAR_tend_NAME(I_cp_rhoq_t(iv)) = trim(AQ_NAME(iv))//'_t_CP'
       VAR_tend_DESC(I_cp_rhoq_t(iv)) = 'tendency rho*'//trim(AQ_NAME(iv))//'in CP'
       VAR_tend_UNIT(I_cp_rhoq_t(iv)) = 'kg/m3/s'
    end do
    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_CP_VARS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_CP_VARS. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_CP_VARS)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** [ATMOS_PHY_CP] prognostic/diagnostic variables'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A15,A,A32,3(A))') &
               '***       |','VARNAME        ','|', 'DESCRIPTION                     ','[', 'UNIT            ',']'
    do iv = 1, VMAX
       if( IO_L ) write(IO_FID_LOG,'(1x,A,i3,A,A15,A,A32,3(A))') &
                  '*** NO.',iv,'|',VAR_NAME(iv),'|',VAR_DESC(iv),'[',VAR_UNIT(iv),']'
    enddo
    ! tendency
    do iv = 1, VMAX_tend
       if( IO_L ) write(IO_FID_LOG,'(1x,A,i3,A,A15,A,A32,3(A))') &
                  '*** NO.',iv+VMAX,'|',VAR_tend_NAME(iv),'|',VAR_tend_DESC(iv),'[',VAR_tend_UNIT(iv),']'
    enddo
    !
    if( IO_L ) write(IO_FID_LOG,*)
    if ( ATMOS_PHY_CP_RESTART_IN_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart input?  : ', trim(ATMOS_PHY_CP_RESTART_IN_BASENAME)
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart input?  : NO'
    endif
    if (       ATMOS_PHY_CP_RESTART_OUTPUT             &
         .AND. ATMOS_PHY_CP_RESTART_OUT_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart output? : ', trim(ATMOS_PHY_CP_RESTART_OUT_BASENAME)
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
    implicit none

    integer :: i,j
    !---------------------------------------------------------------------------
    ! 3d var
    do j  = JS, JE
    do i  = IS, IE
       ATMOS_PHY_CP_cldfrac_dp(   1:KS-1,i,j) = ATMOS_PHY_CP_cldfrac_dp(KS,i,j)
       ATMOS_PHY_CP_cldfrac_dp(KE+1:KA,  i,j) = ATMOS_PHY_CP_cldfrac_dp(KE,i,j)
       ATMOS_PHY_CP_cldfrac_sh(   1:KS-1,i,j) = ATMOS_PHY_CP_cldfrac_sh(KS,i,j)
       ATMOS_PHY_CP_cldfrac_sh(KE+1:KA,  i,j) = ATMOS_PHY_CP_cldfrac_sh(KE,i,j)
       ATMOS_PHY_CP_kf_w0avg(   1:KS-1,i,j)   = ATMOS_PHY_CP_kf_w0avg(KS,i,j)
       ATMOS_PHY_CP_kf_w0avg(KE+1:KA,  i,j)   = ATMOS_PHY_CP_kf_w0avg(KE,i,j)
    enddo
    enddo

    call COMM_vars8( ATMOS_PHY_CP_MFLX_cloudbase (:,:)  , 1 )
    call COMM_vars8( ATMOS_PHY_CP_SFLX_convrain  (:,:)  , 2 )
    call COMM_vars8( ATMOS_PHY_CP_cloudtop       (:,:)  , 3 )
    call COMM_vars8( ATMOS_PHY_CP_cloudbase      (:,:)  , 4 )
    call COMM_vars8( ATMOS_PHY_CP_cldfrac_dp     (:,:,:), 5 )
    call COMM_vars8( ATMOS_PHY_CP_cldfrac_sh     (:,:,:), 6 )
    call COMM_vars8( ATMOS_PHY_CP_kf_nca         (:,:)  , 7 )
    call COMM_vars8( ATMOS_PHY_CP_kf_w0avg       (:,:,:), 8 )

    call COMM_wait ( ATMOS_PHY_CP_MFLX_cloudbase (:,:)  , 1 )
    call COMM_wait ( ATMOS_PHY_CP_SFLX_convrain  (:,:)  , 2 )
    call COMM_wait ( ATMOS_PHY_CP_cloudtop       (:,:)  , 3 )
    call COMM_wait ( ATMOS_PHY_CP_cloudbase      (:,:)  , 4 )
    call COMM_wait ( ATMOS_PHY_CP_cldfrac_dp     (:,:,:), 5 )
    call COMM_wait ( ATMOS_PHY_CP_cldfrac_sh     (:,:,:), 6 )
    call COMM_wait ( ATMOS_PHY_CP_kf_nca         (:,:)  , 7 )
    call COMM_wait ( ATMOS_PHY_CP_kf_w0avg       (:,:,:), 8 )
    ! tendency
    call COMM_vars8( ATMOS_PHY_CP_DENS_t (:,:,:),1+8)
    call COMM_vars8( ATMOS_PHY_CP_RHOT_t (:,:,:),2+8)
    do i = 1,QA
       call COMM_vars8( ATMOS_PHY_CP_RHOQ_t (:,:,:,i),2+i+8)
    end do
    call COMM_wait( ATMOS_PHY_CP_DENS_t (:,:,:),1+8)
    call COMM_wait( ATMOS_PHY_CP_RHOT_t (:,:,:),2+8)
    do i = 1,QA
       call COMM_wait( ATMOS_PHY_CP_RHOQ_t (:,:,:,i),2+i+8)
    end do
    return
  end subroutine ATMOS_PHY_CP_vars_fillhalo

  !-----------------------------------------------------------------------------
  !> Read restart
  subroutine ATMOS_PHY_CP_vars_restart_read
    use scale_fileio, only: &
       FILEIO_read
    use scale_rm_statistics, only: &
       STAT_total
    implicit none
    integer  :: iq
    real(RP) :: total
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input restart file (ATMOS_PHY_CP) ***'

    if ( ATMOS_PHY_CP_RESTART_IN_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(ATMOS_PHY_CP_RESTART_IN_BASENAME)

       call FILEIO_read( ATMOS_PHY_CP_MFLX_cloudbase(:,:), &                            ! [OUT]
                         ATMOS_PHY_CP_RESTART_IN_BASENAME, VAR_NAME(1), 'XY', step=1 )  ! [IN]
       call FILEIO_read( ATMOS_PHY_CP_SFLX_convrain(:,:),  &                            ! [OUT]
                         ATMOS_PHY_CP_RESTART_IN_BASENAME, VAR_NAME(2), 'XY', step=1 )  ! [IN]
       call FILEIO_read( ATMOS_PHY_CP_cloudtop(:,:),       &                            ! [OUT]
                         ATMOS_PHY_CP_RESTART_IN_BASENAME, VAR_NAME(3), 'XY', step=1 )  ! [IN]
       call FILEIO_read( ATMOS_PHY_CP_cloudbase(:,:),      &                            ! [OUT]
                         ATMOS_PHY_CP_RESTART_IN_BASENAME, VAR_NAME(4), 'XY', step=1 )  ! [IN]
       call FILEIO_read( ATMOS_PHY_CP_cldfrac_dp(:,:,:),   &                            ! [OUT]
                         ATMOS_PHY_CP_RESTART_IN_BASENAME, VAR_NAME(5), 'ZXY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_CP_cldfrac_sh(:,:,:),   &                            ! [OUT]
                         ATMOS_PHY_CP_RESTART_IN_BASENAME, VAR_NAME(6), 'ZXY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_CP_kf_nca(:,:),         &                            ! [OUT]
                         ATMOS_PHY_CP_RESTART_IN_BASENAME, VAR_NAME(7), 'XY', step=1 )  ! [IN]
       call FILEIO_read( ATMOS_PHY_CP_kf_w0avg(:,:,:),     &                            ! [OUT]
                         ATMOS_PHY_CP_RESTART_IN_BASENAME, VAR_NAME(8), 'ZXY', step=1 ) ! [IN]
       ! tendency
       call FILEIO_read( ATMOS_PHY_CP_DENS_t(:,:,:),                           &  ! [OUT]
                         ATMOS_PHY_CP_RESTART_IN_BASENAME, VAR_tend_NAME(I_cp_dens_t), 'ZXY', step=1 )  ! [IN]
       call FILEIO_read( ATMOS_PHY_CP_RHOT_t(:,:,:),                           &  ! [OUT]
                         ATMOS_PHY_CP_RESTART_IN_BASENAME, VAR_tend_NAME(I_cp_rhot_t), 'ZXY', step=1 )  ! [IN]
       do iq = 1,QA
          call FILEIO_read( ATMOS_PHY_CP_RHOQ_t(:,:,:,iq),                           &  ! [OUT]
                            ATMOS_PHY_CP_RESTART_IN_BASENAME, VAR_tend_NAME(I_cp_rhoq_t(iq)), 'ZXY', step=1 )  ! [IN]
       end do

       call ATMOS_PHY_CP_vars_fillhalo

       call STAT_total( total, ATMOS_PHY_CP_MFLX_cloudbase(:,:), VAR_NAME(1) )
       call STAT_total( total, ATMOS_PHY_CP_SFLX_convrain(:,:) , VAR_NAME(2) )
       call STAT_total( total, ATMOS_PHY_CP_cloudtop(:,:)      , VAR_NAME(3) )
       call STAT_total( total, ATMOS_PHY_CP_cloudbase(:,:)     , VAR_NAME(4) )
       call STAT_total( total, ATMOS_PHY_CP_cldfrac_dp(:,:,:)  , VAR_NAME(5) )
       call STAT_total( total, ATMOS_PHY_CP_cldfrac_sh(:,:,:)  , VAR_NAME(6) )
       call STAT_total( total, ATMOS_PHY_CP_kf_nca(:,:)        , VAR_NAME(7) )
       call STAT_total( total, ATMOS_PHY_CP_kf_w0avg(:,:,:)    , VAR_NAME(8) )
       ! tendency
       call STAT_total( total, ATMOS_PHY_CP_DENS_t(:,:,:)    , VAR_tend_NAME(I_cp_dens_t) )
       call STAT_total( total, ATMOS_PHY_CP_RHOT_t(:,:,:)    , VAR_tend_NAME(I_cp_rhot_t) )
       do iq = 1,QA
          call STAT_total( total, ATMOS_PHY_CP_RHOQ_t(:,:,:,iq)    , VAR_tend_NAME(I_cp_rhoq_t(iq)) )
       end do
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** restart file for ATMOS_PHY_CP is not specified.'
    endif

    return
  end subroutine ATMOS_PHY_CP_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Write restart
  subroutine ATMOS_PHY_CP_vars_restart_write
    use scale_time, only: &
       TIME_gettimelabel
    use scale_fileio, only: &
       FILEIO_write
    implicit none

    character(len=20)     :: timelabel
    character(len=H_LONG) :: basename
    integer :: iq
    !---------------------------------------------------------------------------

    if ( ATMOS_PHY_CP_RESTART_OUT_BASENAME /= '' ) then

       call TIME_gettimelabel( timelabel )
       write(basename,'(A,A,A)') trim(ATMOS_PHY_CP_RESTART_OUT_BASENAME), '_', trim(timelabel)

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Output restart file (ATMOS_PHY_CP) ***'
       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(basename)

       call FILEIO_write( ATMOS_PHY_CP_MFLX_cloudbase(:,:), basename,  ATMOS_PHY_CP_RESTART_OUT_TITLE, &  ! [IN]
                          VAR_NAME(1), VAR_DESC(1), VAR_UNIT(1), 'XY', ATMOS_PHY_CP_RESTART_OUT_DTYPE  )  ! [IN]
       call FILEIO_write( ATMOS_PHY_CP_SFLX_convrain(:,:), basename,  ATMOS_PHY_CP_RESTART_OUT_TITLE,  &  ! [IN]
                          VAR_NAME(2), VAR_DESC(2), VAR_UNIT(2), 'XY', ATMOS_PHY_CP_RESTART_OUT_DTYPE  )  ! [IN]
       call FILEIO_write( ATMOS_PHY_CP_cloudtop(:,:), basename,  ATMOS_PHY_CP_RESTART_OUT_TITLE,       &  ! [IN]
                          VAR_NAME(3), VAR_DESC(3), VAR_UNIT(3), 'XY', ATMOS_PHY_CP_RESTART_OUT_DTYPE  )  ! [IN]
       call FILEIO_write( ATMOS_PHY_CP_cloudbase(:,:), basename,  ATMOS_PHY_CP_RESTART_OUT_TITLE,      &  ! [IN]
                          VAR_NAME(4), VAR_DESC(4), VAR_UNIT(4), 'XY', ATMOS_PHY_CP_RESTART_OUT_DTYPE  )  ! [IN]
       call FILEIO_write( ATMOS_PHY_CP_cldfrac_dp(:,:,:), basename,  ATMOS_PHY_CP_RESTART_OUT_TITLE,   &  ! [IN]
                          VAR_NAME(5), VAR_DESC(5), VAR_UNIT(5), 'ZXY', ATMOS_PHY_CP_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( ATMOS_PHY_CP_cldfrac_sh(:,:,:), basename,  ATMOS_PHY_CP_RESTART_OUT_TITLE,   &  ! [IN]
                          VAR_NAME(6), VAR_DESC(6), VAR_UNIT(6), 'ZXY', ATMOS_PHY_CP_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( ATMOS_PHY_CP_kf_nca(:,:), basename,  ATMOS_PHY_CP_RESTART_OUT_TITLE,         &  ! [IN]
                          VAR_NAME(7), VAR_DESC(7), VAR_UNIT(7), 'XY', ATMOS_PHY_CP_RESTART_OUT_DTYPE  )  ! [IN]
       call FILEIO_write( ATMOS_PHY_CP_kf_w0avg(:,:,:), basename,  ATMOS_PHY_CP_RESTART_OUT_TITLE,     &  ! [IN]
                          VAR_NAME(8), VAR_DESC(8), VAR_UNIT(8), 'ZXY', ATMOS_PHY_CP_RESTART_OUT_DTYPE  ) ! [IN]
       ! tendency
       call FILEIO_write( ATMOS_PHY_CP_DENS_t(:,:,:), basename,  ATMOS_PHY_CP_RESTART_OUT_TITLE,               & ! [IN]
                       VAR_tend_NAME(I_cp_dens_t), VAR_tend_DESC(I_cp_dens_t), VAR_tend_UNIT(I_cp_dens_t),     & ! [IN]
                      'ZXY', ATMOS_PHY_CP_RESTART_OUT_DTYPE  )                                                   ! [IN]
       call FILEIO_write( ATMOS_PHY_CP_RHOT_t(:,:,:), basename,  ATMOS_PHY_CP_RESTART_OUT_TITLE,               & ! [IN]
                       VAR_tend_NAME(I_cp_rhot_t), VAR_tend_DESC(I_cp_rhot_t), VAR_tend_UNIT(I_cp_rhot_t),     & ! [IN]
                       'ZXY', ATMOS_PHY_CP_RESTART_OUT_DTYPE  )                                                  ! [IN]
       do iq = 1,QA
          call FILEIO_write( ATMOS_PHY_CP_RHOQ_t(:,:,:,iq), basename,  ATMOS_PHY_CP_RESTART_OUT_TITLE,         & ! [IN]
               VAR_tend_NAME(I_cp_rhoq_t(iq)), VAR_tend_DESC(I_cp_rhoq_t(iq)), VAR_tend_UNIT(I_cp_rhoq_t(iq)), & ! [IN]
               'ZXY', ATMOS_PHY_CP_RESTART_OUT_DTYPE  )                                                          ! [IN]
       end do
    endif

    return
  end subroutine ATMOS_PHY_CP_vars_restart_write

  !-----------------------------------------------------------------------------
  !> Create restart file
  subroutine ATMOS_PHY_CP_vars_restart_create
    use scale_time, only: &
       TIME_gettimelabel
    use scale_fileio, only: &
       FILEIO_create
    implicit none

    character(len=20)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    if ( ATMOS_PHY_CP_RESTART_OUT_BASENAME /= '' ) then

       call TIME_gettimelabel( timelabel )
       write(basename,'(A,A,A)') trim(ATMOS_PHY_CP_RESTART_OUT_BASENAME), '_', trim(timelabel)

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Output restart file (ATMOS_PHY_CP) ***'
       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(basename)

       call FILEIO_create(restart_fid,basename,  ATMOS_PHY_CP_RESTART_OUT_TITLE, & ! [IN]
                          ATMOS_PHY_CP_RESTART_OUT_DTYPE  ) ! [IN]

    endif

    return
  end subroutine ATMOS_PHY_CP_vars_restart_create

  !-----------------------------------------------------------------------------
  !> Exit netCDF define mode
  subroutine ATMOS_PHY_CP_vars_restart_enddef
    use scale_fileio, only: &
       FILEIO_enddef
    implicit none

    if ( restart_fid .NE. -1 ) then
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

    if ( restart_fid .NE. -1 ) then
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
    implicit none

    !---------------------------------------------------------------------------

    if ( restart_fid .NE. -1 ) then
       call FILEIO_def_var( restart_fid, VAR_ID(1), VAR_NAME(1), VAR_DESC(1), &
                            VAR_UNIT(1), 'XY', ATMOS_PHY_CP_RESTART_OUT_DTYPE  ) ! [IN]
    endif

    return
  end subroutine ATMOS_PHY_CP_vars_restart_def_var

  !-----------------------------------------------------------------------------
  !> Write restart
  subroutine ATMOS_PHY_CP_vars_restart_write_var
    use scale_fileio, only: &
       FILEIO_write_var
    implicit none

    !---------------------------------------------------------------------------

    if ( restart_fid .NE. -1 ) then
       call FILEIO_write_var( restart_fid, VAR_ID(1), ATMOS_PHY_CP_MFLX_cloudbase(:,:), &
                              VAR_NAME(1), 'XY' ) ! [IN]
    endif

    return
  end subroutine ATMOS_PHY_CP_vars_restart_write_var

end module mod_atmos_phy_cp_vars
