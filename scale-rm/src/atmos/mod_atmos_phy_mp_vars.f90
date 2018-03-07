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
!! @li      2015-09-08 (Y.Sato)       [add] Add ATMOS_PHY_MP_EVAPORATE
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
  use scale_atmos_grid_cartesC_index
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

  public :: ATMOS_PHY_MP_vars_restart_create
  public :: ATMOS_PHY_MP_vars_restart_open
  public :: ATMOS_PHY_MP_vars_restart_def_var
  public :: ATMOS_PHY_MP_vars_restart_enddef
  public :: ATMOS_PHY_MP_vars_restart_close

  public :: ATMOS_PHY_MP_vars_history

  public :: ATMOS_PHY_MP_vars_get_diagnostic
  public :: ATMOS_PHY_MP_vars_reset_diagnostics

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical,                public :: ATMOS_PHY_MP_RESTART_OUTPUT                = .false.                !< output restart file?

  character(len=H_LONG),  public :: ATMOS_PHY_MP_RESTART_IN_BASENAME           = ''                     !< Basename of the input  file
  logical,                public :: ATMOS_PHY_MP_RESTART_IN_AGGREGATE                                   !< Switch to use aggregate file
  logical,                public :: ATMOS_PHY_MP_RESTART_IN_POSTFIX_TIMELABEL  = .false.                !< Add timelabel to the basename of input  file?
  character(len=H_LONG),  public :: ATMOS_PHY_MP_RESTART_OUT_BASENAME          = ''                     !< Basename of the output file
  logical,                public :: ATMOS_PHY_MP_RESTART_OUT_AGGREGATE                                  !< Switch to use aggregate file
  logical,                public :: ATMOS_PHY_MP_RESTART_OUT_POSTFIX_TIMELABEL = .true.                 !< Add timelabel to the basename of output file?
  character(len=H_MID),   public :: ATMOS_PHY_MP_RESTART_OUT_TITLE             = 'ATMOS_PHY_MP restart' !< title    of the output file
  character(len=H_SHORT), public :: ATMOS_PHY_MP_RESTART_OUT_DTYPE             = 'DEFAULT'              !< REAL4 or REAL8

  real(RP), public :: ATMOS_PHY_MP_cldfrac_thleshold

  real(RP), public, allocatable :: ATMOS_PHY_MP_DENS_t(:,:,:)    ! tendency DENS [kg/m3/s]
  real(RP), public, allocatable :: ATMOS_PHY_MP_MOMZ_t(:,:,:)    ! tendency MOMZ [kg/m2/s2]
  real(RP), public, allocatable :: ATMOS_PHY_MP_RHOU_t(:,:,:)    ! tendency dens*U [kg/m2/s2]
  real(RP), public, allocatable :: ATMOS_PHY_MP_RHOV_t(:,:,:)    ! tendency dens*V [kg/m2/s2]
  real(RP), public, allocatable :: ATMOS_PHY_MP_RHOT_t(:,:,:)    ! tendency RHOT [K*kg/m3/s]
  real(RP), public, allocatable :: ATMOS_PHY_MP_RHOQ_t(:,:,:,:)  ! tendency rho*QTRC [kg/kg/s]
  real(RP), public, allocatable :: ATMOS_PHY_MP_RHOH  (:,:,:)    ! diabatic heating rate [J/kg/s]

  ! obsolute
  real(RP), public, allocatable :: ATMOS_PHY_MP_MOMX_t(:,:,:)    ! tendency MOMX [kg/m2/s2]
  real(RP), public, allocatable :: ATMOS_PHY_MP_MOMY_t(:,:,:)    ! tendency MOMY [kg/m2/s2]

  real(RP), public, allocatable :: ATMOS_PHY_MP_EVAPORATE(:,:,:) ! number concentration of evaporated cloud [/m3]
  real(RP), public, allocatable :: ATMOS_PHY_MP_SFLX_rain(:,:)   ! precipitation flux (liquid) [kg/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_MP_SFLX_snow(:,:)   ! precipitation flux (solid)  [kg/m2/s]

  integer, public :: QA_MP
  integer, public :: QS_MP
  integer, public :: QE_MP

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
  integer,                private            :: VAR_ID(VMAX)   !< ID    of the variables
  integer,                private            :: restart_fid = -1  ! file ID

  data VAR_NAME / 'SFLX_rain', &
                  'SFLX_snow'  /
  data VAR_DESC / 'precipitation flux (liquid)', &
                  'precipitation flux (solid)'   /
  data VAR_UNIT / 'kg/m2/s', &
                  'kg/m2/s' /


  ! for diagnostics
  real(RP), private, allocatable :: ATMOS_PHY_MP_CLDFRAC(:,:,:)
  real(RP), private, allocatable :: ATMOS_PHY_MP_Re     (:,:,:,:)
  real(RP), private, allocatable :: ATMOS_PHY_MP_Qe     (:,:,:,:)
  logical, private :: DIAG_CLDFRAC
  logical, private :: DIAG_Re
  logical, private :: DIAG_Qe


  ! for history
  integer, private             :: HIST_CLDFRAC_id
  integer, private,allocatable :: HIST_Re_id(:)
  logical, private             :: HIST_Re


  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_MP_vars_setup
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_atmos_hydrometeor, only: &
       N_HYD, &
       QHA,  &
       HYD_NAME, &
       HYD_DESC
    use scale_file_history, only: &
       FILE_HISTORY_reg
    implicit none

    NAMELIST / PARAM_ATMOS_PHY_MP_VARS / &
       ATMOS_PHY_MP_RESTART_IN_BASENAME,           &
       ATMOS_PHY_MP_RESTART_IN_AGGREGATE,          &
       ATMOS_PHY_MP_RESTART_IN_POSTFIX_TIMELABEL,  &
       ATMOS_PHY_MP_RESTART_OUTPUT,                &
       ATMOS_PHY_MP_RESTART_OUT_BASENAME,          &
       ATMOS_PHY_MP_RESTART_OUT_AGGREGATE,         &
       ATMOS_PHY_MP_RESTART_OUT_POSTFIX_TIMELABEL, &
       ATMOS_PHY_MP_RESTART_OUT_TITLE,             &
       ATMOS_PHY_MP_RESTART_OUT_DTYPE

    integer :: ierr
    integer :: iv, ih
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[VARS] / Categ[ATMOS PHY_MP] / Origin[SCALE-RM]'

    allocate( ATMOS_PHY_MP_DENS_t   (KA,IA,JA)    )
    allocate( ATMOS_PHY_MP_MOMZ_t   (KA,IA,JA)    )
    allocate( ATMOS_PHY_MP_RHOU_t   (KA,IA,JA)    )
    allocate( ATMOS_PHY_MP_RHOV_t   (KA,IA,JA)    )
    allocate( ATMOS_PHY_MP_RHOT_t   (KA,IA,JA)    )
    allocate( ATMOS_PHY_MP_RHOQ_t   (KA,IA,JA,QS_MP:QE_MP) )
    allocate( ATMOS_PHY_MP_RHOH     (KA,IA,JA)    )
    allocate( ATMOS_PHY_MP_EVAPORATE(KA,IA,JA)    )
    allocate( ATMOS_PHY_MP_MOMX_t   (KA,IA,JA)    )
    allocate( ATMOS_PHY_MP_MOMY_t   (KA,IA,JA)    )
    ! tentative approach
    ATMOS_PHY_MP_DENS_t   (:,:,:)   = 0.0_RP
    ATMOS_PHY_MP_MOMZ_t   (:,:,:)   = 0.0_RP
    ATMOS_PHY_MP_RHOU_t   (:,:,:)   = 0.0_RP
    ATMOS_PHY_MP_RHOV_t   (:,:,:)   = 0.0_RP
    ATMOS_PHY_MP_RHOT_t   (:,:,:)   = 0.0_RP
    ATMOS_PHY_MP_RHOQ_t   (:,:,:,:) = 0.0_RP
    ATMOS_PHY_MP_RHOH     (:,:,:)   = 0.0_RP
    ATMOS_PHY_MP_EVAPORATE(:,:,:)   = 0.0_RP
    ATMOS_PHY_MP_MOMX_t   (:,:,:)   = 0.0_RP
    ATMOS_PHY_MP_MOMY_t   (:,:,:)   = 0.0_RP

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
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_ATMOS_PHY_MP_VARS)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** [ATMOS_PHY_MP] prognostic/diagnostic variables'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A24,A,A48,A,A12,A)') &
               '***       |', 'VARNAME                 ','|', &
               'DESCRIPTION                                     ', '[', 'UNIT        ', ']'
    do iv = 1, VMAX
       if( IO_L ) write(IO_FID_LOG,'(1x,A,I3,A,A24,A,A48,A,A12,A)') &
                  '*** NO.',iv,'|',VAR_NAME(iv),'|',VAR_DESC(iv),'[',VAR_UNIT(iv),']'
    enddo

    if( IO_L ) write(IO_FID_LOG,*)
    if ( ATMOS_PHY_MP_RESTART_IN_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart input?  : YES, file = ', trim(ATMOS_PHY_MP_RESTART_IN_BASENAME)
       if( IO_L ) write(IO_FID_LOG,*) '*** Add timelabel?  : ', ATMOS_PHY_MP_RESTART_IN_POSTFIX_TIMELABEL
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart input?  : NO'
    endif
    if (       ATMOS_PHY_MP_RESTART_OUTPUT             &
         .AND. ATMOS_PHY_MP_RESTART_OUT_BASENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart output? : YES, file = ', trim(ATMOS_PHY_MP_RESTART_OUT_BASENAME)
       if( IO_L ) write(IO_FID_LOG,*) '*** Add timelabel?  : ', ATMOS_PHY_MP_RESTART_OUT_POSTFIX_TIMELABEL
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Restart output? : NO'
       ATMOS_PHY_MP_RESTART_OUTPUT = .false.
    endif


    ! diagnostices
    allocate( ATMOS_PHY_MP_CLDFRAC(KA,IA,JA) )
    allocate( ATMOS_PHY_MP_Re     (KA,IA,JA,N_HYD) )
    allocate( ATMOS_PHY_MP_Qe     (KA,IA,JA,N_HYD) )
!OCL XFILL
    ATMOS_PHY_MP_CLDFRAC(:,:,:) = UNDEF
!OCL XFILL
    ATMOS_PHY_MP_Re     (:,:,:,:) = UNDEF
!OCL XFILL
    ATMOS_PHY_MP_Qe     (:,:,:,:) = UNDEF
    DIAG_CLDFRAC = .false.
    DIAG_Re      = .false.
    DIAG_Qe      = .false.

    ! history
    call FILE_HISTORY_reg( 'CLDFRAC', 'cloud fraction', '1', HIST_CLDFRAC_id, fill_halo=.true., dim_type='ZXY' )

    HIST_Re = .false.
    allocate( HIST_Re_id(N_HYD) )
    do ih = 1, N_HYD
       call FILE_HISTORY_reg( 'Re_'//trim(HYD_NAME(ih)), 'effective radius of '//trim(HYD_DESC(ih)), 'cm', HIST_Re_id(ih), fill_halo=.true., dim_type='ZXY' )
       if ( HIST_Re_id(ih) > 0 ) HIST_Re = .true.
    end do

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
  !> Open restart file for read
  subroutine ATMOS_PHY_MP_vars_restart_open
    use scale_time, only: &
       TIME_gettimelabel
    use scale_file_cartesC, only: &
       FILE_CARTESC_open
    implicit none

    character(len=19)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Open restart file (ATMOS_PHY_MP) ***'

    if ( ATMOS_PHY_MP_RESTART_IN_BASENAME /= '' ) then

       if ( ATMOS_PHY_MP_RESTART_IN_POSTFIX_TIMELABEL ) then
          call TIME_gettimelabel( timelabel )
          basename = trim(ATMOS_PHY_MP_RESTART_IN_BASENAME)//'_'//trim(timelabel)
       else
          basename = trim(ATMOS_PHY_MP_RESTART_IN_BASENAME)
       endif

       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(basename)

       call FILE_CARTESC_open( basename, restart_fid, aggregate=ATMOS_PHY_MP_RESTART_IN_AGGREGATE )
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** restart file for ATMOS_PHY_MP is not specified.'
    endif

    return
  end subroutine ATMOS_PHY_MP_vars_restart_open

  !-----------------------------------------------------------------------------
  !> Read restart
  subroutine ATMOS_PHY_MP_vars_restart_read
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STATISTICS_total
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_AREA, &
       ATMOS_GRID_CARTESC_REAL_TOTAREA
    use scale_file, only: &
       FILE_get_aggregate
    use scale_file_cartesC, only: &
       FILE_CARTESC_read, &
       FILE_CARTESC_flush
    implicit none

    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Read from restart file (ATMOS_PHY_MP) ***'

       call FILE_CARTESC_read( restart_fid, VAR_NAME(1), 'XY', & ! [IN]
                               ATMOS_PHY_MP_SFLX_rain(:,:)     ) ! [OUT]

       call FILE_CARTESC_read( restart_fid, VAR_NAME(2), 'XY', & ! [IN]
                               ATMOS_PHY_MP_SFLX_snow(:,:)     ) ! [OUT]

       if ( FILE_get_AGGREGATE(restart_fid) ) then
          call FILE_CARTESC_flush( restart_fid ) ! X/Y halos have been read from file
       else
          call ATMOS_PHY_MP_vars_fillhalo
       end if

       if ( STATISTICS_checktotal ) then
          call STATISTICS_total( IA, IS, IE, JA, JS, JE, &
                                 ATMOS_PHY_MP_SFLX_rain(:,:), VAR_NAME(1), & ! (in)
                                 ATMOS_GRID_CARTESC_REAL_AREA(:,:),        & ! (in)
                                 ATMOS_GRID_CARTESC_REAL_TOTAREA           ) ! (in)
          call STATISTICS_total( IA, IS, IE, JA, JS, JE, &
                                 ATMOS_PHY_MP_SFLX_snow(:,:), VAR_NAME(2), & ! (in)
                                 ATMOS_GRID_CARTESC_REAL_AREA(:,:),        & ! (in)
                                 ATMOS_GRID_CARTESC_REAL_TOTAREA           ) ! (in)
       endif
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** invalid restart file ID for ATMOS_PHY_MP.'
    endif

    return
  end subroutine ATMOS_PHY_MP_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Create restart file
  subroutine ATMOS_PHY_MP_vars_restart_create
    use scale_time, only: &
       TIME_gettimelabel
    use scale_file_cartesC, only: &
       FILE_CARTESC_create
    implicit none

    character(len=19)     :: timelabel
    character(len=H_LONG) :: basename
    !---------------------------------------------------------------------------

    if ( ATMOS_PHY_MP_RESTART_OUT_BASENAME /= '' ) then

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Create restart file (ATMOS_PHY_AE) ***'

       if ( ATMOS_PHY_MP_RESTART_OUT_POSTFIX_TIMELABEL ) then
          call TIME_gettimelabel( timelabel )
          basename = trim(ATMOS_PHY_MP_RESTART_OUT_BASENAME)//'_'//trim(timelabel)
       else
          basename = trim(ATMOS_PHY_MP_RESTART_OUT_BASENAME)
       endif

       if( IO_L ) write(IO_FID_LOG,*) '*** basename: ', trim(basename)

       call FILE_CARTESC_create( &
            basename, ATMOS_PHY_MP_RESTART_OUT_TITLE, ATMOS_PHY_MP_RESTART_OUT_DTYPE, & ! [IN]
            restart_fid,                                                              & ! [OUT]
            aggregate=ATMOS_PHY_MP_RESTART_OUT_AGGREGATE                              ) ! [IN]

    endif

    return
  end subroutine ATMOS_PHY_MP_vars_restart_create

  !-----------------------------------------------------------------------------
  !> Exit netCDF define mode
  subroutine ATMOS_PHY_MP_vars_restart_enddef
    use scale_file_cartesC, only: &
       FILE_CARTESC_enddef
    implicit none

    if ( restart_fid /= -1 ) then
       call FILE_CARTESC_enddef( restart_fid ) ! [IN]
    endif

    return
  end subroutine ATMOS_PHY_MP_vars_restart_enddef

  !-----------------------------------------------------------------------------
  !> Close restart file
  subroutine ATMOS_PHY_MP_vars_restart_close
    use scale_file_cartesC, only: &
       FILE_CARTESC_close
    implicit none
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Close restart file (ATMOS_PHY_MP) ***'

       call FILE_CARTESC_close( restart_fid ) ! [IN]

       restart_fid = -1
    endif

    return
  end subroutine ATMOS_PHY_MP_vars_restart_close

  !-----------------------------------------------------------------------------
  !> Define variables in restart file
  subroutine ATMOS_PHY_MP_vars_restart_def_var
    use scale_file_cartesC, only: &
       FILE_CARTESC_def_var
    implicit none
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then

       call FILE_CARTESC_def_var( restart_fid, VAR_NAME(1), VAR_DESC(1), VAR_UNIT(1), 'XY', ATMOS_PHY_MP_RESTART_OUT_DTYPE, &
                                  VAR_ID(1) )
       call FILE_CARTESC_def_var( restart_fid, VAR_NAME(2), VAR_DESC(2), VAR_UNIT(2), 'XY', ATMOS_PHY_MP_RESTART_OUT_DTYPE, &
                                  VAR_ID(2) )
    endif

    return
  end subroutine ATMOS_PHY_MP_vars_restart_def_var

  !-----------------------------------------------------------------------------
  !> Write restart
  subroutine ATMOS_PHY_MP_vars_restart_write
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STATISTICS_total
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_AREA, &
       ATMOS_GRID_CARTESC_REAL_TOTAREA
    use scale_file_cartesC, only: &
       FILE_CARTESC_write_var
    implicit none
    !---------------------------------------------------------------------------

    if ( restart_fid /= -1 ) then

       call ATMOS_PHY_MP_vars_fillhalo

       if ( STATISTICS_checktotal ) then
          call STATISTICS_total( IA, IS, IE, JA, JS, JE, &
                                 ATMOS_PHY_MP_SFLX_rain(:,:), VAR_NAME(1), & ! (in)
                                 ATMOS_GRID_CARTESC_REAL_AREA(:,:),        & ! (in)
                                 ATMOS_GRID_CARTESC_REAL_TOTAREA           ) ! (in)
          call STATISTICS_total( IA, IS, IE, JA, JS, JE, &
                                 ATMOS_PHY_MP_SFLX_snow(:,:), VAR_NAME(2), & ! (in)
                                 ATMOS_GRID_CARTESC_REAL_AREA(:,:),        & ! (in)
                                 ATMOS_GRID_CARTESC_REAL_TOTAREA           ) ! (in)
       endif

       call FILE_CARTESC_write_var( restart_fid, VAR_ID(1), ATMOS_PHY_MP_SFLX_rain(:,:), &
                              VAR_NAME(1), 'XY' ) ! [IN]
       call FILE_CARTESC_write_var( restart_fid, VAR_ID(2), ATMOS_PHY_MP_SFLX_snow(:,:), &
                              VAR_NAME(2), 'XY' ) ! [IN]

    endif

    return
  end subroutine ATMOS_PHY_MP_vars_restart_write

  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_vars_history( &
       DENS, TEMP, QTRC )
    use scale_atmos_hydrometeor, only: &
       N_HYD
    use scale_file_history, only: &
       FILE_HISTORY_query, &
       FILE_HISTORY_put
    implicit none

    real(RP), intent(in) :: DENS(KA,IA,JA)
    real(RP), intent(in) :: TEMP(KA,IA,JA)
    real(RP), intent(in) :: QTRC(KA,IA,JA,QA)

    real(RP) :: WORK  (KA,IA,JA,N_HYD)
    logical  :: do_put
    integer  :: ih
    !---------------------------------------------------------------------------

    if ( HIST_CLDFRAC_id > 0 ) then
       call FILE_HISTORY_query( HIST_CLDFRAC_id, do_put )

       if ( do_put ) then
          call ATMOS_PHY_MP_vars_get_diagnostic( &
               DENS(:,:,:), TEMP(:,:,:), QTRC(:,:,:,:), & ! [IN]
               CLDFRAC=WORK(:,:,:,1)                    ) ! [OUT]
          call FILE_HISTORY_put( HIST_CLDFRAC_id, WORK(:,:,:,1) )
       end if
    end if

    if ( HIST_Re ) then
       do ih = 1, N_HYD
          if ( HIST_Re_id(ih) > 0 ) then
             call FILE_HISTORY_query( HIST_Re_id(ih), do_put )
             if ( do_put ) then
                call ATMOS_PHY_MP_vars_get_diagnostic( &
                     DENS(:,:,:), TEMP(:,:,:), QTRC(:,:,:,:), & ! [IN]
                     Re=WORK(:,:,:,:)                         ) ! [OUT]
                exit
             end if
          end if
       end do
       if ( do_put ) then
          do ih = 1, N_HYD
             if ( HIST_Re_id(ih) > 0 ) then
                call FILE_HISTORY_query( HIST_Re_id(ih), do_put )
                if ( do_put ) call FILE_HISTORY_put( HIST_Re_id(ih), WORK(:,:,:,ih) )
             end if
          end do
       end if
    end if

    return
  end subroutine ATMOS_PHY_MP_vars_history

  subroutine ATMOS_PHY_MP_vars_get_diagnostic( &
       DENS, TEMP, QTRC, &
       CLDFRAC, Re, Qe   )
    use scale_atmos_hydrometeor, only: &
       N_HYD, &
       I_HC,  &
       I_HR,  &
       I_HI,  &
       I_HS,  &
       I_HG,  &
       I_HH,  &
       QHS,   &
       QHE
    use scale_atmos_phy_mp_kessler, only: &
       ATMOS_PHY_MP_KESSLER_mass_ratio, &
       ATMOS_PHY_MP_KESSLER_effective_radius, &
       ATMOS_PHY_MP_KESSLER_cloud_fraction
    use scale_atmos_phy_mp_tomita08, only: &
       ATMOS_PHY_MP_TOMITA08_mass_ratio, &
       ATMOS_PHY_MP_TOMITA08_effective_radius, &
       ATMOS_PHY_MP_TOMITA08_cloud_fraction
    use scale_atmos_phy_mp, only: &
       ATMOS_PHY_MP_CloudFraction,   &
       ATMOS_PHY_MP_EffectiveRadius, &
       ATMOS_PHY_MP_MixingRatio
    use mod_atmos_admin, only: &
       ATMOS_PHY_MP_TYPE

    real(RP), intent(in) :: DENS(KA,IA,JA)
    real(RP), intent(in) :: TEMP(KA,IA,JA)
    real(RP), intent(in) :: QTRC(KA,IA,JA,QA)
    real(RP), intent(out), optional :: CLDFRAC(KA,IA,JA)       !> cloud fraction [0-1]
    real(RP), intent(out), optional :: Re     (KA,IA,JA,N_HYD) !> effective radius [cm]
    real(RP), intent(out), optional :: Qe     (KA,IA,JA,N_HYD) !> mass ratio [kg/kg]

    integer :: k, i, j, ih

    if ( present(CLDFRAC) ) then
       if ( .not. DIAG_CLDFRAC ) then
          select case ( ATMOS_PHY_MP_TYPE )
          case ( 'KESSLER' )
             call ATMOS_PHY_MP_kessler_cloud_fraction( &
                  KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                  QTRC(:,:,:,QHS:QHE), ATMOS_PHY_MP_cldfrac_thleshold, & ! [IN]
                  ATMOS_PHY_MP_CLDFRAC(:,:,:)                          ) ! [OUT]
          case ( 'TOMITA08' )
             call ATMOS_PHY_MP_tomita08_cloud_fraction( &
                  KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                  QTRC(:,:,:,QHS:QHE), ATMOS_PHY_MP_cldfrac_thleshold, & ! [IN]
                  ATMOS_PHY_MP_CLDFRAC(:,:,:)                          ) ! [OUT]
          case default
             if ( associated(ATMOS_PHY_MP_CloudFraction) ) then
                call ATMOS_PHY_MP_CloudFraction( &
                     ATMOS_PHY_MP_CLDFRAC(:,:,:),                   & ! [OUT]
                     QTRC(:,:,:,:), ATMOS_PHY_MP_cldfrac_thleshold  ) ! [IN]
             else
!OCL XFILL
                ATMOS_PHY_MP_CLDFRAC(:,:,:) = 0.0_RP
             end if
          end select
          DIAG_CLDFRAC = .true.
       end if
!OCL XFILL
       do j = JSB, JEB
       do i = ISB, IEB
       do k = KS, KE
          CLDFRAC(k,i,j) = ATMOS_PHY_MP_CLDFRAC(k,i,j)
       end do
       end do
       end do
    end if

    if ( present(Re) ) then
       if ( .not. DIAG_Re ) then
          select case ( ATMOS_PHY_MP_TYPE )
          case ( 'KESSLER' )
             call ATMOS_PHY_MP_kessler_effective_radius( &
                  KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                  DENS(:,:,:), TEMP(:,:,:), QTRC(:,:,:,QHS:QHE), & ! [IN]
                  ATMOS_PHY_MP_Re(:,:,:,:)                       ) ! [OUT]
          case ( 'TOMITA08' )
             call ATMOS_PHY_MP_tomita08_effective_radius( &
                  KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                  DENS(:,:,:), TEMP(:,:,:), QTRC(:,:,:,QHS:QHE), & ! [IN]
                  ATMOS_PHY_MP_Re(:,:,:,:)                       ) ! [OUT]
          case default
             if ( associated(ATMOS_PHY_MP_EffectiveRadius) ) then
                call ATMOS_PHY_MP_EffectiveRadius( &
                     ATMOS_PHY_MP_Re(:,:,:,:),               & ! [OUT]
                     QTRC(:,:,:,:), DENS(:,:,:), TEMP(:,:,:) ) ! [IN]
             else
!OCL XFILL
                ATMOS_PHY_MP_Re(:,:,:,:) = 0.0_RP
             end if
          end select
          DIAG_Re = .true.
       end if
!OCL XFILL
       do ih = 1, N_HYD
       do j = JSB, JEB
       do i = ISB, IEB
       do k = KS, KE
          Re(k,i,j,ih) = ATMOS_PHY_MP_Re(k,i,j,ih)
       end do
       end do
       end do
       end do
    end if

    if ( present(Qe) ) then
       if ( .not. DIAG_Qe ) then
          select case ( ATMOS_PHY_MP_TYPE )
          case ( 'KESSLER' )
             call ATMOS_PHY_MP_kessler_mass_ratio( &
                  KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                  QTRC(:,:,:,QHS:QHE),     & ! [IN]
                  ATMOS_PHY_MP_Qe(:,:,:,:) ) ! [OUT]
          case ( 'TOMITA08' )
             call ATMOS_PHY_MP_tomita08_mass_ratio( &
                  KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                  QTRC(:,:,:,QHS:QHE),     & ! [IN]
                  ATMOS_PHY_MP_Qe(:,:,:,:) ) ! [OUT]
          case default
             if ( associated(ATMOS_PHY_MP_MixingRatio) ) then
                call ATMOS_PHY_MP_MixingRatio( &
                     ATMOS_PHY_MP_Qe(:,:,:,:),  & ! [OUT]
                     QTRC(:,:,:,:)              ) ! [IN]
             else
!OCL XIFLL
                ATMOS_PHY_MP_Qe(:,:,:,:) = 0.0_RP
             end if
          end select
          DIAG_Qe = .true.
       end if
!OCL XIFLL
       do ih = 1, N_HYD
       do j = JSB, JEB
       do i = ISB, IEB
       do k = KS, KE
          Qe(k,i,j,ih) = ATMOS_PHY_MP_Qe(k,i,j,ih)
       end do
       end do
       end do
       end do
    end if

    return
  end subroutine ATMOS_PHY_MP_vars_get_diagnostic

  subroutine ATMOS_PHY_MP_vars_reset_diagnostics
    DIAG_CLDFRAC = .false.
    DIAG_Re      = .false.
    DIAG_Qe      = .false.

    return
  end subroutine ATMOS_PHY_MP_vars_reset_diagnostics

end module mod_atmos_phy_mp_vars
