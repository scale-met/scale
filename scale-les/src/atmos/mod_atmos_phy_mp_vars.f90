!-------------------------------------------------------------------------------
!> module Atmosphere / Physics Cloud Microphysics
!!
!! @par Description
!!          Container for mod_atmos_phy_mp
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-05-04 (H.Yashiro)    [new]
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module mod_atmos_phy_mp_vars
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
  public :: ATMOS_PHY_MP_vars_setup
  public :: ATMOS_PHY_MP_vars_fillhalo
  public :: ATMOS_PHY_MP_vars_restart_read
  public :: ATMOS_PHY_MP_vars_restart_write

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical,               public :: ATMOS_PHY_MP_RESTART_OUTPUT       = .false.                !< output restart file?

  character(len=H_LONG), public :: ATMOS_PHY_MP_RESTART_IN_BASENAME  = ''                     !< basename of the restart file
  character(len=H_LONG), public :: ATMOS_PHY_MP_RESTART_OUT_BASENAME = ''                     !< basename of the output file
  character(len=H_MID),  public :: ATMOS_PHY_MP_RESTART_OUT_TITLE    = 'ATMOS_PHY_MP restart' !< title    of the output file
  character(len=H_MID),  public :: ATMOS_PHY_MP_RESTART_OUT_DTYPE    = 'DEFAULT'              !< REAL4 or REAL8

  real(RP), public, allocatable :: ATMOS_PHY_MP_DENS_t(:,:,:)   ! tendency DENS [kg/m3/s]
  real(RP), public, allocatable :: ATMOS_PHY_MP_MOMZ_t(:,:,:)   ! tendency MOMZ [kg/m2/s2]
  real(RP), public, allocatable :: ATMOS_PHY_MP_MOMX_t(:,:,:)   ! tendency MOMX [kg/m2/s2]
  real(RP), public, allocatable :: ATMOS_PHY_MP_MOMY_t(:,:,:)   ! tendency MOMY [kg/m2/s2]
  real(RP), public, allocatable :: ATMOS_PHY_MP_RHOT_t(:,:,:)   ! tendency RHOT [K*kg/m3/s]
  real(RP), public, allocatable :: ATMOS_PHY_MP_RHOQ_t(:,:,:,:) ! tendency rho*QTRC [kg/kg/s]

  real(RP), public, allocatable :: ATMOS_PHY_MP_SFLX_rain(:,:) ! precipitation flux (liquid) [kg/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_MP_SFLX_snow(:,:) ! precipitation flux (solid)  [kg/m2/s]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,                private, parameter :: VMAX = 2       !< number of the variables
  integer,                private, parameter :: I_SFLX_rain = 1
  integer,                private, parameter :: I_SFLX_snow = 2

  character(len=H_SHORT), private            :: VAR_NAME(VMAX) !< name  of the variables
  character(len=H_MID),   private            :: VAR_DESC(VMAX) !< desc. of the variables
  character(len=H_SHORT), private            :: VAR_UNIT(VMAX) !< unit  of the variables

  data VAR_NAME / 'SFLX_rain', &
                  'SFLX_snow'  /
  data VAR_DESC / 'precipitation flux (liquid)', &
                  'precipitation flux (solid)'   /
  data VAR_UNIT / 'kg/m2/s', &
                  'kg/m2/s' /

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_MP_vars_setup
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    implicit none

    NAMELIST / PARAM_ATMOS_PHY_MP_VARS / &
       ATMOS_PHY_MP_RESTART_IN_BASENAME,  &
       ATMOS_PHY_MP_RESTART_OUTPUT,       &
       ATMOS_PHY_MP_RESTART_OUT_BASENAME, &
       ATMOS_PHY_MP_RESTART_OUT_TITLE,    &
       ATMOS_PHY_MP_RESTART_OUT_DTYPE

    integer :: ierr
    integer :: iv
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[VARS] / Categ[ATMOS_PHY_MP] / Origin[SCALE-LES]'

    allocate( ATMOS_PHY_MP_DENS_t(KA,IA,JA)    )
    allocate( ATMOS_PHY_MP_MOMZ_t(KA,IA,JA)    )
    allocate( ATMOS_PHY_MP_MOMX_t(KA,IA,JA)    )
    allocate( ATMOS_PHY_MP_MOMY_t(KA,IA,JA)    )
    allocate( ATMOS_PHY_MP_RHOT_t(KA,IA,JA)    )
    allocate( ATMOS_PHY_MP_RHOQ_t(KA,IA,JA,QA) )
    ! tentative approach
    ATMOS_PHY_MP_DENS_t(:,:,:)   = 0.0_RP
    ATMOS_PHY_MP_MOMZ_t(:,:,:)   = 0.0_RP
    ATMOS_PHY_MP_MOMX_t(:,:,:)   = 0.0_RP
    ATMOS_PHY_MP_MOMY_t(:,:,:)   = 0.0_RP
    ATMOS_PHY_MP_RHOT_t(:,:,:)   = 0.0_RP
    ATMOS_PHY_MP_RHOQ_t(:,:,:,:) = 0.0_RP

    allocate( ATMOS_PHY_MP_SFLX_rain(IA,JA) )
    allocate( ATMOS_PHY_MP_SFLX_snow(IA,JA) )
    ATMOS_PHY_MP_SFLX_rain(:,:) = UNDEF
    ATMOS_PHY_MP_SFLX_snow(:,:) = UNDEF

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_MP_VARS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_MP_VARS. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_MP_VARS)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** [ATMOS_PHY_MP] prognostic/diagnostic variables'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A15,A,A32,3(A))') &
               '***       |','VARNAME        ','|', 'DESCRIPTION                     ','[', 'UNIT            ',']'
    do iv = 1, VMAX
       if( IO_L ) write(IO_FID_LOG,'(1x,A,i3,A,A15,A,A32,3(A))') &
                  '*** NO.',iv,'|',VAR_NAME(iv),'|',VAR_DESC(iv),'[',VAR_UNIT(iv),']'
    enddo

    if( IO_L ) write(IO_FID_LOG,*)
    if ( ATMOS_PHY_MP_RESTART_IN_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart input?  : ', trim(ATMOS_PHY_MP_RESTART_IN_BASENAME)
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart input?  : NO'
    endif
    if (       ATMOS_PHY_MP_RESTART_OUTPUT             &
         .AND. ATMOS_PHY_MP_RESTART_OUT_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart output? : ', trim(ATMOS_PHY_MP_RESTART_OUT_BASENAME)
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart output? : NO'
       ATMOS_PHY_MP_RESTART_OUTPUT = .false.
    endif

    return
  end subroutine ATMOS_PHY_MP_vars_setup

  !-----------------------------------------------------------------------------
  !> HALO Communication
  subroutine ATMOS_PHY_MP_vars_fillhalo
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none
    !---------------------------------------------------------------------------

    call COMM_vars8( ATMOS_PHY_MP_SFLX_rain(:,:), 1 )
    call COMM_vars8( ATMOS_PHY_MP_SFLX_snow(:,:), 2 )
    call COMM_wait ( ATMOS_PHY_MP_SFLX_rain(:,:), 1 )
    call COMM_wait ( ATMOS_PHY_MP_SFLX_snow(:,:), 2 )

    return
  end subroutine ATMOS_PHY_MP_vars_fillhalo

  !-----------------------------------------------------------------------------
  !> Read restart
  subroutine ATMOS_PHY_MP_vars_restart_read
    use scale_fileio, only: &
       FILEIO_read
    use scale_statistics, only: &
       STAT_total
    implicit none

    real(RP) :: total
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input restart file (ATMOS_PHY_MP) ***'

    if ( ATMOS_PHY_MP_RESTART_IN_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(ATMOS_PHY_MP_RESTART_IN_BASENAME)

       call FILEIO_read( ATMOS_PHY_MP_SFLX_rain(:,:),                                & ! [OUT]
                         ATMOS_PHY_MP_RESTART_IN_BASENAME, VAR_NAME(1), 'XY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_MP_SFLX_snow(:,:),                                & ! [OUT]
                         ATMOS_PHY_MP_RESTART_IN_BASENAME, VAR_NAME(2), 'XY', step=1 ) ! [IN]

       call ATMOS_PHY_MP_vars_fillhalo

       call STAT_total( total, ATMOS_PHY_MP_SFLX_rain(:,:), VAR_NAME(1) )
       call STAT_total( total, ATMOS_PHY_MP_SFLX_snow(:,:), VAR_NAME(2) )
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** restart file for ATMOS_PHY_MP is not specified.'
    endif

    return
  end subroutine ATMOS_PHY_MP_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Write restart
  subroutine ATMOS_PHY_MP_vars_restart_write
    use scale_time, only: &
       TIME_gettimelabel
    use scale_fileio, only: &
       FILEIO_write
    use scale_statistics, only: &
       STAT_total
    implicit none

    character(len=15)     :: timelabel
    character(len=H_LONG) :: basename

    real(RP) :: total
    !---------------------------------------------------------------------------

    if ( ATMOS_PHY_MP_RESTART_OUT_BASENAME /= '' ) then

       call TIME_gettimelabel( timelabel )
       write(basename,'(A,A,A)') trim(ATMOS_PHY_MP_RESTART_OUT_BASENAME), '_', trim(timelabel)

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Output restart file (ATMOS_PHY_MP) ***'
       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(basename)

       call ATMOS_PHY_MP_vars_fillhalo

       call STAT_total( total, ATMOS_PHY_MP_SFLX_rain(:,:), VAR_NAME(1) )
       call STAT_total( total, ATMOS_PHY_MP_SFLX_snow(:,:), VAR_NAME(1) )

       call FILEIO_write( ATMOS_PHY_MP_SFLX_rain(:,:), basename,       ATMOS_PHY_MP_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(1), VAR_DESC(1), VAR_UNIT(1), 'XY', ATMOS_PHY_MP_RESTART_OUT_DTYPE  ) ! [IN]
       call FILEIO_write( ATMOS_PHY_MP_SFLX_snow(:,:), basename,       ATMOS_PHY_MP_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(2), VAR_DESC(2), VAR_UNIT(2), 'XY', ATMOS_PHY_MP_RESTART_OUT_DTYPE  ) ! [IN]

    endif

    return
  end subroutine ATMOS_PHY_MP_vars_restart_write

end module mod_atmos_phy_mp_vars
