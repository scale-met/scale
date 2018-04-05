!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Aerosol Microphysics
!!
!! @par Description
!!          Container for mod_atmos_phy_ae
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_atmos_phy_ae_vars
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_atmos_grid_cartesC_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_AE_vars_setup
  public :: ATMOS_PHY_AE_vars_fillhalo
  public :: ATMOS_PHY_AE_vars_restart_read
  public :: ATMOS_PHY_AE_vars_restart_write

  public :: ATMOS_PHY_AE_vars_restart_create
  public :: ATMOS_PHY_AE_vars_restart_open
  public :: ATMOS_PHY_AE_vars_restart_def_var
  public :: ATMOS_PHY_AE_vars_restart_enddef
  public :: ATMOS_PHY_AE_vars_restart_close

  public :: ATMOS_PHY_AE_vars_history

  public :: ATMOS_PHY_AE_vars_get_diagnostic
  public :: ATMOS_PHY_AE_vars_reset_diagnostics

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical,                public :: ATMOS_PHY_AE_RESTART_OUTPUT                = .false.                !< output restart file?

  character(len=H_LONG),  public :: ATMOS_PHY_AE_RESTART_IN_BASENAME           = ''                     !< Basename of the input  file
  logical,                public :: ATMOS_PHY_AE_RESTART_IN_AGGREGATE                                   !< Switch to use aggregate file
  logical,                public :: ATMOS_PHY_AE_RESTART_IN_POSTFIX_TIMELABEL  = .false.                !< Add timelabel to the basename of input  file?
  character(len=H_LONG),  public :: ATMOS_PHY_AE_RESTART_OUT_BASENAME          = ''                     !< Basename of the output file
  logical,                public :: ATMOS_PHY_AE_RESTART_OUT_AGGREGATE                                   !< Switch to use aggregate file
  logical,                public :: ATMOS_PHY_AE_RESTART_OUT_POSTFIX_TIMELABEL = .true.                 !< Add timelabel to the basename of output file?
  character(len=H_MID),   public :: ATMOS_PHY_AE_RESTART_OUT_TITLE             = 'ATMOS_PHY_AE restart' !< title    of the output file
  character(len=H_SHORT), public :: ATMOS_PHY_AE_RESTART_OUT_DTYPE             = 'DEFAULT'              !< REAL4 or REAL8

  real(RP), public, allocatable :: ATMOS_PHY_AE_RHOQ_t(:,:,:,:)                               ! tendency rho*QTRC [kg/kg/s]

  real(RP), public, allocatable :: ATMOS_PHY_AE_CCN(:,:,:)                                    ! cloud condensation nuclei [/m3]
  real(RP), public, allocatable :: ATMOS_PHY_AE_CCN_t(:,:,:)                                  ! tendency CCN [/m3/s]
  real(RP), public, allocatable :: ATMOS_PHY_AE_EMIT(:,:,:,:)                                 ! emission of aerosol and gas

  integer, public :: QA_AE
  integer, public :: QS_AE
  integer, public :: QE_AE

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,                private, parameter :: VMAX  = 1       !< number of the variables
  integer,                private, parameter :: I_CCN = 1

  character(len=H_SHORT), private            :: VAR_NAME(VMAX) !< name  of the variables
  character(len=H_MID),   private            :: VAR_DESC(VMAX) !< desc. of the variables
  character(len=H_SHORT), private            :: VAR_UNIT(VMAX) !< unit  of the variables
  integer,                private            :: VAR_ID(VMAX)   !< ID    of the variables
  integer,                private            :: restart_fid = -1  ! file ID

  data VAR_NAME / 'CCN' /
  data VAR_DESC / 'cloud condensation nuclei' /
  data VAR_UNIT / 'num/m3' /


  ! for diagnostics
  real(RP), private, allocatable :: ATMOS_PHY_AE_Re(:,:,:,:)
  real(RP), private, allocatable :: ATMOS_PHY_AE_Qe(:,:,:,:)
  logical, private :: DIAG_Re
  logical, private :: DIAG_Qe

  ! for history
  integer, private,allocatable :: HIST_Re_id(:)
  logical, private             :: HIST_Re

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_AE_vars_setup
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_atmos_aerosol, only: &
       N_AE, &
       AE_NAME, &
       AE_DESC
    use scale_file_history, only: &
       FILE_HISTORY_reg
    implicit none

    NAMELIST / PARAM_ATMOS_PHY_AE_VARS / &
       ATMOS_PHY_AE_RESTART_IN_BASENAME,           &
       ATMOS_PHY_AE_RESTART_IN_AGGREGATE,          &
       ATMOS_PHY_AE_RESTART_IN_POSTFIX_TIMELABEL,  &
       ATMOS_PHY_AE_RESTART_OUTPUT,                &
       ATMOS_PHY_AE_RESTART_OUT_BASENAME,          &
       ATMOS_PHY_AE_RESTART_OUT_AGGREGATE,         &
       ATMOS_PHY_AE_RESTART_OUT_POSTFIX_TIMELABEL, &
       ATMOS_PHY_AE_RESTART_OUT_TITLE,             &
       ATMOS_PHY_AE_RESTART_OUT_DTYPE

    integer :: ierr
    integer :: iv
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_PROGRESS(*) 'Module[VARS] / Categ[ATMOS PHY_AE] / Origin[SCALE-RM]'

    allocate( ATMOS_PHY_AE_RHOQ_t(KA,IA,JA,QS_AE:QE_AE) )
    ATMOS_PHY_AE_RHOQ_t(:,:,:, :) = UNDEF

    allocate( ATMOS_PHY_AE_CCN(KA,IA,JA) )
    ATMOS_PHY_AE_CCN(:,:,:) = UNDEF

    allocate( ATMOS_PHY_AE_CCN_t(KA,IA,JA) )
    ATMOS_PHY_AE_CCN_t(:,:,:) = UNDEF

    allocate( ATMOS_PHY_AE_EMIT(KA,IA,JA,QS_AE:QE_AE) )
    ATMOS_PHY_AE_EMIT(:,:,:,:) = 0.0_RP

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_AE_VARS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_PHY_AE_vars_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_AE_vars_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_AE_VARS. Check!'
       call PRC_abort
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_ATMOS_PHY_AE_VARS)

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_AE_vars_setup",*) '[ATMOS_PHY_AE] prognostic/diagnostic variables'
    LOG_INFO_CONT('(1x,A,A24,A,A48,A,A12,A)') &
               '***       |', 'VARNAME                 ','|', &
               'DESCRIPTION                                     ', '[', 'UNIT        ', ']'
    do iv = 1, VMAX
       LOG_INFO_CONT('(1x,A,I3,A,A24,A,A48,A,A12,A)') &
                  '*** NO.',iv,'|',VAR_NAME(iv),'|',VAR_DESC(iv),'[',VAR_UNIT(iv),']'
    enddo

    LOG_NEWLINE
    if ( ATMOS_PHY_AE_RESTART_IN_BASENAME /= '' ) then
       LOG_INFO("ATMOS_PHY_AE_vars_setup",*) 'Restart input?  : YES, file = ', trim(ATMOS_PHY_AE_RESTART_IN_BASENAME)
       LOG_INFO("ATMOS_PHY_AE_vars_setup",*) 'Add timelabel?  : ', ATMOS_PHY_AE_RESTART_IN_POSTFIX_TIMELABEL
    else
       LOG_INFO("ATMOS_PHY_AE_vars_setup",*) 'Restart input?  : NO'
    endif
    if (       ATMOS_PHY_AE_RESTART_OUTPUT             &
         .AND. ATMOS_PHY_AE_RESTART_OUT_BASENAME /= '' ) then
       LOG_INFO("ATMOS_PHY_AE_vars_setup",*) 'Restart output? : YES, file = ', trim(ATMOS_PHY_AE_RESTART_OUT_BASENAME)
       LOG_INFO("ATMOS_PHY_AE_vars_setup",*) 'Add timelabel?  : ', ATMOS_PHY_AE_RESTART_OUT_POSTFIX_TIMELABEL
    else
       LOG_INFO("ATMOS_PHY_AE_vars_setup",*) 'Restart output? : NO'
       ATMOS_PHY_AE_RESTART_OUTPUT = .false.
    endif

    ! diagnostices
    allocate( ATMOS_PHY_AE_Re(KA,IA,JA,N_AE) )
    allocate( ATMOS_PHY_AE_Qe(KA,IA,JA,N_AE) )
!OCL XFILL
    ATMOS_PHY_AE_Re(:,:,:,:) = UNDEF
!OCL XFILL
    ATMOS_PHY_AE_Qe(:,:,:,:) = UNDEF
    DIAG_Re = .false.
    DIAG_Qe = .false.

    ! history
    HIST_Re = .false.
    allocate( HIST_Re_id(N_AE) )
    do iv = 1, N_AE
       call FILE_HISTORY_reg( 'Re_'//trim(AE_NAME(iv)), 'effective radius of '//trim(AE_DESC(iv)), 'cm', HIST_Re_id(iv), fill_halo=.true., dim_type='ZXY' )
       if ( HIST_Re_id(iv) > 0 ) HIST_Re = .true.
    end do

    return
  end subroutine ATMOS_PHY_AE_vars_setup

  !-----------------------------------------------------------------------------
  !> HALO Communication
  subroutine ATMOS_PHY_AE_vars_fillhalo
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    integer :: i, j
    !---------------------------------------------------------------------------

    do j  = JS, JE
    do i  = IS, IE
       ATMOS_PHY_AE_CCN(   1:KS-1,i,j) = ATMOS_PHY_AE_CCN(KS,i,j)
       ATMOS_PHY_AE_CCN(KE+1:KA,  i,j) = ATMOS_PHY_AE_CCN(KE,i,j)
    enddo
    enddo

    call COMM_vars8( ATMOS_PHY_AE_CCN(:,:,:), 1 )
    call COMM_wait ( ATMOS_PHY_AE_CCN(:,:,:), 1 )

    return
  end subroutine ATMOS_PHY_AE_vars_fillhalo

  !-----------------------------------------------------------------------------
  !> Open restart file for read
  subroutine ATMOS_PHY_AE_vars_restart_open
    use scale_time, only: &
       TIME_gettimelabel
    use scale_file_cartesC, only: &
       FILE_CARTESC_open
    implicit none

    character(len=19)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_AE_vars_restart_open",*) 'Open restart file (ATMOS_PHY_AE) '

    if ( ATMOS_PHY_AE_RESTART_IN_BASENAME /= '' ) then

       if ( ATMOS_PHY_AE_RESTART_IN_POSTFIX_TIMELABEL ) then
          call TIME_gettimelabel( timelabel )
          basename = trim(ATMOS_PHY_AE_RESTART_IN_BASENAME)//'_'//trim(timelabel)
       else
          basename = trim(ATMOS_PHY_AE_RESTART_IN_BASENAME)
       endif

       LOG_INFO("ATMOS_PHY_AE_vars_restart_open",*) 'basename: ', trim(basename)

       call FILE_CARTESC_open( basename, restart_fid, aggregate=ATMOS_PHY_AE_RESTART_IN_AGGREGATE )
    else
       LOG_INFO("ATMOS_PHY_AE_vars_restart_open",*) 'restart file for ATMOS_PHY_AE is not specified.'
    endif

    return
  end subroutine ATMOS_PHY_AE_vars_restart_open

  !-----------------------------------------------------------------------------
  !> Read restart
  subroutine ATMOS_PHY_AE_vars_restart_read
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STATISTICS_total
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_VOL, &
       ATMOS_GRID_CARTESC_REAL_TOTVOL
    use scale_file, only: &
       FILE_get_aggregate
    use scale_file_cartesC, only: &
       FILE_CARTESC_read, &
       FILE_CARTESC_flush
    implicit none

    integer  :: i, j
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then
       LOG_NEWLINE
       LOG_INFO("ATMOS_PHY_AE_vars_restart_read",*) 'Read from restart file (ATMOS_PHY_AE) '

       call FILE_CARTESC_read( restart_fid, VAR_NAME(1), 'ZXY', & ! [IN]
                               ATMOS_PHY_AE_CCN(:,:,:)          ) ! [OUT]
                         

       if ( FILE_get_AGGREGATE(restart_fid) ) then
          call FILE_CARTESC_flush( restart_fid ) ! X/Y halos have been read from file

          ! fill K halos
          do j  = 1, JA
          do i  = 1, IA
             ATMOS_PHY_AE_CCN(   1:KS-1,i,j) = ATMOS_PHY_AE_CCN(KS,i,j)
             ATMOS_PHY_AE_CCN(KE+1:KA,  i,j) = ATMOS_PHY_AE_CCN(KE,i,j)
          enddo
          enddo
       else
          call ATMOS_PHY_AE_vars_fillhalo
       end if

       if ( STATISTICS_checktotal ) then
          call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                 ATMOS_PHY_AE_CCN(:,:,:), VAR_NAME(1), &
                                 ATMOS_GRID_CARTESC_REAL_VOL(:,:,:),   & ! (in)
                                 ATMOS_GRID_CARTESC_REAL_TOTVOL        ) ! (in)
       end if
    else
       LOG_INFO("ATMOS_PHY_AE_vars_restart_read",*) 'invlaid restart file ID for ATMOS_PHY_AE.'
    endif

    return
  end subroutine ATMOS_PHY_AE_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Create restart file
  subroutine ATMOS_PHY_AE_vars_restart_create
    use scale_time, only: &
       TIME_gettimelabel
    use scale_file_cartesC, only: &
       FILE_CARTESC_create
    implicit none

    character(len=19)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    if ( ATMOS_PHY_AE_RESTART_OUT_BASENAME /= '' ) then

       LOG_NEWLINE
       LOG_INFO("ATMOS_PHY_AE_vars_restart_create",*) 'Create restart file (ATMOS_PHY_AE) '

       if ( ATMOS_PHY_AE_RESTART_OUT_POSTFIX_TIMELABEL ) then
          call TIME_gettimelabel( timelabel )
          basename = trim(ATMOS_PHY_AE_RESTART_OUT_BASENAME)//'_'//trim(timelabel)
       else
          basename = trim(ATMOS_PHY_AE_RESTART_OUT_BASENAME)
       endif

       LOG_INFO("ATMOS_PHY_AE_vars_restart_create",*) 'basename: ', trim(basename)

       call FILE_CARTESC_create( &
            basename, ATMOS_PHY_AE_RESTART_OUT_TITLE, ATMOS_PHY_AE_RESTART_OUT_DTYPE, & ! [IN]
            restart_fid,                                                              & ! [OUT]
            aggregate=ATMOS_PHY_AE_RESTART_OUT_AGGREGATE                              ) ! [IN]

    endif

    return
  end subroutine ATMOS_PHY_AE_vars_restart_create

  !-----------------------------------------------------------------------------
  !> Exit netCDF define mode
  subroutine ATMOS_PHY_AE_vars_restart_enddef
    use scale_file_cartesC, only: &
       FILE_CARTESC_enddef
    implicit none

    if ( restart_fid /= -1 ) then
       call FILE_CARTESC_enddef( restart_fid ) ! [IN]
    endif

    return
  end subroutine ATMOS_PHY_AE_vars_restart_enddef

  !-----------------------------------------------------------------------------
  !> Close restart file
  subroutine ATMOS_PHY_AE_vars_restart_close
    use scale_file_cartesC, only: &
       FILE_CARTESC_close
    implicit none
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then
       LOG_NEWLINE
       LOG_INFO("ATMOS_PHY_AE_vars_restart_close",*) 'Close restart file (ATMOS_PHY_AE) '

       call FILE_CARTESC_close( restart_fid ) ! [IN]

       restart_fid = -1
    endif

    return
  end subroutine ATMOS_PHY_AE_vars_restart_close

  !-----------------------------------------------------------------------------
  !> Write restart
  subroutine ATMOS_PHY_AE_vars_restart_def_var
    use scale_file_cartesC, only: &
       FILE_CARTESC_def_var
    implicit none
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then
       call FILE_CARTESC_def_var( restart_fid, VAR_NAME(1), VAR_DESC(1), VAR_UNIT(1), 'ZXY', ATMOS_PHY_AE_RESTART_OUT_DTYPE, &
                                  VAR_ID(1) )
    endif

    return
  end subroutine ATMOS_PHY_AE_vars_restart_def_var

  !-----------------------------------------------------------------------------
  !> Write restart
  subroutine ATMOS_PHY_AE_vars_restart_write
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STATISTICS_total
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_VOL, &
       ATMOS_GRID_CARTESC_REAL_TOTVOL
    use scale_file_cartesC, only: &
       FILE_CARTESC_write_var
    implicit none
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then

       call ATMOS_PHY_AE_vars_fillhalo

       if ( STATISTICS_checktotal ) then
          call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                 ATMOS_PHY_AE_CCN(:,:,:), VAR_NAME(1), &
                                 ATMOS_GRID_CARTESC_REAL_VOL(:,:,:),   & ! (in)
                                 ATMOS_GRID_CARTESC_REAL_TOTVOL        ) ! (in)
       end if

       call FILE_CARTESC_write_var( restart_fid, VAR_ID(1), ATMOS_PHY_AE_CCN(:,:,:), VAR_NAME(1), 'ZXY' ) ! [IN]

    endif

    return
  end subroutine ATMOS_PHY_AE_vars_restart_write

  subroutine ATMOS_PHY_AE_vars_history( &
       QTRC, RH )
    use scale_tracer, only: &
       QA
    use scale_atmos_aerosol, only: &
       N_AE
    use scale_file_history, only: &
       FILE_HISTORY_put
    real(RP), intent(in) :: QTRC(KA,IA,JA,QA)
    real(RP), intent(in) :: RH(KA,IA,JA)

    real(RP) :: WORK(KA,IA,JA,N_AE)
    integer  :: iv

    if ( HIST_Re ) then
       call ATMOS_PHY_AE_vars_get_diagnostic( &
            QTRC(:,:,:,:), RH(:,:,:), & ! [IN]
            Re=WORK(:,:,:,:)          ) ! [OUT]
       do iv = 1, N_AE
          if ( HIST_Re_id(iv) > 0 ) &
               call FILE_HISTORY_put( HIST_Re_id(iv), WORK(:,:,:,iv) )
       end do
    end if

    return
  end subroutine ATMOS_PHY_AE_vars_history

  subroutine ATMOS_PHY_AE_vars_get_diagnostic( &
       QTRC, RH, &
       Re, Qe   )
    use scale_tracer, only: &
       QA
    use scale_atmos_aerosol, only: &
       N_AE
    use scale_atmos_phy_ae_kajino13, only: &
       ATMOS_PHY_AE_kajino13_effective_radius
    use mod_atmos_admin, only: &
       ATMOS_PHY_AE_TYPE

    real(RP), intent(in) :: QTRC(KA,IA,JA,QA)
    real(RP), intent(in) :: RH(KA,IA,JA)
    real(RP), intent(out), optional :: Re(KA,IA,JA,N_AE) !> effective radius [cm]
    real(RP), intent(out), optional :: Qe(KA,IA,JA,N_AE) !> mass ratio [kg/kg]

    if ( present(Re) ) then
       if ( .not. DIAG_Re ) then
          select case ( ATMOS_PHY_AE_TYPE )
          case ( 'KAJINO13' )
             call ATMOS_PHY_AE_kajino13_effective_radius( &
                     KA, IA, JA, QA_AE, &
                     QTRC(:,:,:,QS_AE:QE_AE), RH(:,:,:), & ! [IN]
                     ATMOS_PHY_AE_Re(:,:,:,:)            ) ! [OUT]
          case default
             ATMOS_PHY_AE_Re(:,:,:,:) = 0.0_RP
          end select
          DIAG_Re = .true.
       end if
!OCL XFILL
       Re(:,:,:,:) = ATMOS_PHY_AE_Re(:,:,:,:)
    end if

    if ( present(Qe) ) then
       if ( .not. DIAG_Qe ) then
          select case ( ATMOS_PHY_AE_TYPE )
          case default
             ATMOS_PHY_AE_Qe(:,:,:,:) = 0.0_RP
          end select
          DIAG_Qe = .true.
       end if
!OCL XIFLL
       Qe(:,:,:,:) = ATMOS_PHY_AE_Qe(:,:,:,:)
    end if

    return
  end subroutine ATMOS_PHY_AE_vars_get_diagnostic

  subroutine ATMOS_PHY_AE_vars_reset_diagnostics
    DIAG_Re      = .false.
    DIAG_Qe      = .false.

    return
  end subroutine ATMOS_PHY_AE_vars_reset_diagnostics

end module mod_atmos_phy_ae_vars
