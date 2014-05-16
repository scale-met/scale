!-------------------------------------------------------------------------------
!> module Atmosphere / Physics Turbulence
!!
!! @par Description
!!          Container for mod_atmos_phy_tb
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-05-04 (H.Yashiro)    [new]
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module mod_atmos_phy_tb_vars
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
  public :: ATMOS_PHY_TB_vars_setup
  public :: ATMOS_PHY_TB_vars_fillhalo
  public :: ATMOS_PHY_TB_vars_restart_read
  public :: ATMOS_PHY_TB_vars_restart_write

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, allocatable :: ATMOS_PHY_TB_MOMZ_t(:,:,:)   ! tendency MOMZ [kg/m2/s2]
  real(RP), public, allocatable :: ATMOS_PHY_TB_MOMX_t(:,:,:)   ! tendency MOMX [kg/m2/s2]
  real(RP), public, allocatable :: ATMOS_PHY_TB_MOMY_t(:,:,:)   ! tendency MOMY [kg/m2/s2]
  real(RP), public, allocatable :: ATMOS_PHY_TB_RHOT_t(:,:,:)   ! tendency RHOT [K*kg/m3/s]
  real(RP), public, allocatable :: ATMOS_PHY_TB_QTRC_t(:,:,:,:) ! tendency QTRC [kg/kg/s]

  real(RP), public, allocatable :: ATMOS_PHY_TB_TKE(:,:,:) ! turburent kinetic energy [m2/s2]
  real(RP), public, allocatable :: ATMOS_PHY_TB_nu (:,:,:) ! eddy viscosity           [m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_TB_Ri (:,:,:) ! Richardson number        [NIL]
  real(RP), public, allocatable :: ATMOS_PHY_TB_Pr (:,:,:) ! Prantle number           [NIL]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical,                private, save      :: ATMOS_PHY_TB_RESTART_OUTPUT        = .false.
  character(len=H_LONG),  private, save      :: ATMOS_PHY_TB_RESTART_IN_BASENAME   = ''
  character(len=H_LONG),  private, save      :: ATMOS_PHY_TB_RESTART_OUT_BASENAME  = ''
  character(len=H_MID),   private, save      :: ATMOS_PHY_TB_RESTART_OUT_TITLE     = 'ATMOS_PHY_TB restart'
  character(len=H_MID),   private, save      :: ATMOS_PHY_TB_RESTART_OUT_DTYPE     = 'DEFAULT'              !< REAL4 or REAL8

  integer,                private, parameter :: VMAX = 4       !< number of the variables
  character(len=H_SHORT), private, save      :: VAR_NAME(VMAX) !< name  of the variables
  character(len=H_MID),   private, save      :: VAR_DESC(VMAX) !< desc. of the variables
  character(len=H_SHORT), private, save      :: VAR_UNIT(VMAX) !< unit  of the variables

  data VAR_NAME / 'TKE', &
                  'nu',  &
                  'Ri',  &
                  'Pr'   /
  data VAR_DESC / 'turburent kinetic energy', &
                  'eddy viscosity',           &
                  'Prantle number',           &
                  'Richardson number'         /
  data VAR_UNIT / 'm2/s2', &
                  'm2/s',  &
                  'NIL',   &
                  'NIL'    /

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_TB_vars_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    NAMELIST / PARAM_ATMOS_PHY_TB_VARS / &
       ATMOS_PHY_TB_RESTART_IN_BASENAME, &
       ATMOS_PHY_TB_RESTART_OUTPUT,      &
       ATMOS_PHY_TB_RESTART_OUT_BASENAME

    integer :: v, ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[ATMOS_PHY_TB VARS]/Categ[ATMOS]'
    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_TB_VARS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_TB_VARS. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_TB_VARS)

    allocate( ATMOS_PHY_TB_MOMZ_t(KA,IA,JA)    )
    allocate( ATMOS_PHY_TB_MOMX_t(KA,IA,JA)    )
    allocate( ATMOS_PHY_TB_MOMY_t(KA,IA,JA)    )
    allocate( ATMOS_PHY_TB_RHOT_t(KA,IA,JA)    )
    allocate( ATMOS_PHY_TB_QTRC_t(KA,IA,JA,QA) )

    allocate( ATMOS_PHY_TB_TKE(KA,IA,JA) )
    allocate( ATMOS_PHY_TB_nu (KA,IA,JA) )
    allocate( ATMOS_PHY_TB_Ri (KA,IA,JA) )
    allocate( ATMOS_PHY_TB_Pr (KA,IA,JA) )

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** [ATMOS_PHY_TB] prognostic/diagnostic variables'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A8,A,A32,3(A))') &
               '***       |',' VARNAME','|', 'DESCRIPTION                     ','[', 'UNIT            ',']'
    do v = 1, VMAX
       if( IO_L ) write(IO_FID_LOG,'(1x,A,i3,A,A8,A,A32,3(A))') &
                  '*** NO.',v,'|',trim(VAR_NAME(v)),'|',VAR_DESC(v),'[',VAR_UNIT(v),']'
    enddo

    ! restart switch
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) 'Output...'
    if ( ATMOS_PHY_TB_RESTART_OUTPUT ) then
       if( IO_L ) write(IO_FID_LOG,*) '  Restart output (ATMOS_PHY_TB) : YES'
    else
       if( IO_L ) write(IO_FID_LOG,*) '  Restart output (ATMOS_PHY_TB) : NO'
    endif
    if( IO_L ) write(IO_FID_LOG,*)

    return
  end subroutine ATMOS_PHY_TB_vars_setup

  !-----------------------------------------------------------------------------
  !> HALO Communication
  subroutine ATMOS_PHY_TB_vars_fillhalo
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none
    !---------------------------------------------------------------------------

    do j  = JS, JE
    do i  = IS, IE
       ATMOS_PHY_TB_TKE(   1:KS-1,i,j) = ATMOS_PHY_TB_TKE(KS,i,j)
       ATMOS_PHY_TB_nu (   1:KS-1,i,j) = ATMOS_PHY_TB_nu (KS,i,j)
       ATMOS_PHY_TB_Ri (   1:KS-1,i,j) = ATMOS_PHY_TB_Ri (KS,i,j)
       ATMOS_PHY_TB_Pr (   1:KS-1,i,j) = ATMOS_PHY_TB_Pr (KS,i,j)
       ATMOS_PHY_TB_TKE(KE+1:KA,  i,j) = ATMOS_PHY_TB_TKE(KE,i,j)
       ATMOS_PHY_TB_nu (KE+1:KA,  i,j) = ATMOS_PHY_TB_nu (KE,i,j)
       ATMOS_PHY_TB_Ri (KE+1:KA,  i,j) = ATMOS_PHY_TB_Ri (KE,i,j)
       ATMOS_PHY_TB_Pr (KE+1:KA,  i,j) = ATMOS_PHY_TB_Pr (KE,i,j)
    enddo
    enddo

    call COMM_vars8( ATMOS_PHY_TB_TKE(:,:,:), 1 )
    call COMM_vars8( ATMOS_PHY_TB_nu (:,:,:), 2 )
    call COMM_vars8( ATMOS_PHY_TB_Ri (:,:,:), 3 )
    call COMM_vars8( ATMOS_PHY_TB_Pr (:,:,:), 4 )
    call COMM_wait ( ATMOS_PHY_TB_TKE(:,:,:), 1 )
    call COMM_wait ( ATMOS_PHY_TB_nu (:,:,:), 2 )
    call COMM_wait ( ATMOS_PHY_TB_Ri (:,:,:), 3 )
    call COMM_wait ( ATMOS_PHY_TB_Pr (:,:,:), 4 )

    return
  end subroutine ATMOS_PHY_TB_vars_fillhalo

  !-----------------------------------------------------------------------------
  !> Read restart
  subroutine ATMOS_PHY_TB_vars_restart_read
    use scale_fileio, only: &
       FILEIO_read
    use scale_const, only: &
       CONST_UNDEF
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input restart file (ATMOS_PHY_TB) ***'

    if ( ATMOS_PHY_TB_RESTART_IN_BASENAME /= '' ) then

       call FILEIO_read( ATMOS_PHY_TB_TKE(:,:,:),                                     & ! [OUT]
                         ATMOS_PHY_TB_RESTART_IN_BASENAME, VAR_NAME(1), 'ZXY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_TB_nu (:,:,:),                                     & ! [OUT]
                         ATMOS_PHY_TB_RESTART_IN_BASENAME, VAR_NAME(2), 'ZXY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_TB_Ri (:,:,:),                                     & ! [OUT]
                         ATMOS_PHY_TB_RESTART_IN_BASENAME, VAR_NAME(3), 'ZXY', step=1 ) ! [IN]
       call FILEIO_read( ATMOS_PHY_TB_Pr (:,:,:),                                     & ! [OUT]
                         ATMOS_PHY_TB_RESTART_IN_BASENAME, VAR_NAME(4), 'ZXY', step=1 ) ! [IN]

       call ATMOS_PHY_TB_vars_fillhalo

    else

       if( IO_L ) write(IO_FID_LOG,*) '*** restart file for ATMOS_PHY_TB is not specified.'

       ATMOS_PHY_TB_TKE(:,:,:) = CONST_UNDEF
       ATMOS_PHY_TB_nu (:,:,:) = CONST_UNDEF
       ATMOS_PHY_TB_Ri (:,:,:) = CONST_UNDEF
       ATMOS_PHY_TB_Pr (:,:,:) = CONST_UNDEF

    endif

    return
  end subroutine ATMOS_PHY_TB_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Write restart
  subroutine ATMOS_PHY_TB_vars_restart_write
    use scale_time, only: &
       TIME_gettimelabel
    use scale_fileio, only: &
       FILEIO_write
    implicit none

    character(len=15)     :: timelabel
    character(len=H_LONG) :: basename

    integer :: n
    !---------------------------------------------------------------------------

    if ( ATMOS_PHY_TB_RESTART_OUT_BASENAME /= '' ) then

       call TIME_gettimelabel( timelabel )
       write(basename,'(A,A,A)') trim(ATMOS_PHY_TB_RESTART_OUT_BASENAME), '_', trim(timelabel)

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Output restart file (ATMOS_PHY_TB) ***'
       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(basename)

       call FILEIO_write( ATMOS_PHY_TB_TKE(:,:,:), basename,            ATMOS_PHY_TB_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(1), VAR_DESC(1), VAR_UNIT(1), 'ZXY', ATMOS_PHY_TB_RESTART_OUT_DTYPE, ) ! [IN]
       call FILEIO_write( ATMOS_PHY_TB_nu (:,:,:), basename,            ATMOS_PHY_TB_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(2), VAR_DESC(2), VAR_UNIT(2), 'ZXY', ATMOS_PHY_TB_RESTART_OUT_DTYPE, ) ! [IN]
       call FILEIO_write( ATMOS_PHY_TB_Ri (:,:,:), basename,            ATMOS_PHY_TB_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(3), VAR_DESC(3), VAR_UNIT(3), 'ZXY', ATMOS_PHY_TB_RESTART_OUT_DTYPE, ) ! [IN]
       call FILEIO_write( ATMOS_PHY_TB_Pr (:,:,:), basename,            ATMOS_PHY_TB_RESTART_OUT_TITLE, & ! [IN]
                          VAR_NAME(4), VAR_DESC(4), VAR_UNIT(4), 'ZXY', ATMOS_PHY_TB_RESTART_OUT_DTYPE, ) ! [IN]

    endif

    return
  end subroutine ATMOS_PHY_TB_vars_restart_write

end module mod_atmos_phy_tb_vars
