!-------------------------------------------------------------------------------
!> module Atmosphere / Physics Cloud Microphysics
!!
!! @par Description
!!          Container for mod_atmos_phy_mp
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_atmos_phy_mp_vars
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_icoA_index
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_MP_vars_setup
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

  real(RP), public, allocatable :: ATMOS_PHY_MP_EVAPORATE(:,:,:,:) ! number concentration of evaporated cloud [/m3]
  real(RP), public, allocatable :: ATMOS_PHY_MP_SFLX_rain(:,:,:)   ! precipitation flux (liquid) [kg/m2/s]
  real(RP), public, allocatable :: ATMOS_PHY_MP_SFLX_snow(:,:,:)   ! precipitation flux (solid)  [kg/m2/s]

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
  real(RP), private, allocatable :: ATMOS_PHY_MP_CLDFRAC(:,:,:,:)
  real(RP), private, allocatable :: ATMOS_PHY_MP_Re     (:,:,:,:,:)
  real(RP), private, allocatable :: ATMOS_PHY_MP_Qe     (:,:,:,:,:)
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
    use scale_prc, only: &
       PRC_abort
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

    namelist / PARAM_ATMOS_PHY_MP_VARS / &
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

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_MP_vars_setup",*) 'Setup'

    allocate( ATMOS_PHY_MP_EVAPORATE(KA,IA,JA,ADM_lall)    )
    ! tentative approach
    ATMOS_PHY_MP_EVAPORATE(:,:,:,:)   = 0.0_RP

    allocate( ATMOS_PHY_MP_SFLX_rain(IA,JA,ADM_lall) )
    allocate( ATMOS_PHY_MP_SFLX_snow(IA,JA,ADM_lall) )
    ATMOS_PHY_MP_SFLX_rain(:,:,:) = UNDEF
    ATMOS_PHY_MP_SFLX_snow(:,:,:) = UNDEF

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_MP_VARS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_PHY_MP_vars_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_MP_vars_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_MP_VARS. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_MP_VARS)

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_MP_vars_setup",*) '[ATMOS_PHY_MP] prognostic/diagnostic variables'
    LOG_INFO_CONT('(1x,A,A24,A,A48,A,A12,A)') &
               '      |', 'VARNAME                 ','|', &
               'DESCRIPTION                                     ', '[', 'UNIT        ', ']'
    do iv = 1, VMAX
       LOG_INFO_CONT('(1x,A,I3,A,A24,A,A48,A,A12,A)') &
                  'NO.',iv,'|',VAR_NAME(iv),'|',VAR_DESC(iv),'[',VAR_UNIT(iv),']'
    enddo

    LOG_NEWLINE
    if ( ATMOS_PHY_MP_RESTART_IN_BASENAME /= '' ) then
       LOG_INFO("ATMOS_PHY_MP_vars_setup",*) 'Restart input?  : YES, file = ', trim(ATMOS_PHY_MP_RESTART_IN_BASENAME)
       LOG_INFO("ATMOS_PHY_MP_vars_setup",*) 'Add timelabel?  : ', ATMOS_PHY_MP_RESTART_IN_POSTFIX_TIMELABEL
    else
       LOG_INFO("ATMOS_PHY_MP_vars_setup",*) 'Restart input?  : NO'
    endif
    if (       ATMOS_PHY_MP_RESTART_OUTPUT             &
         .AND. ATMOS_PHY_MP_RESTART_OUT_BASENAME /= '' ) then
       LOG_INFO("ATMOS_PHY_MP_vars_setup",*) 'Restart output? : YES, file = ', trim(ATMOS_PHY_MP_RESTART_OUT_BASENAME)
       LOG_INFO("ATMOS_PHY_MP_vars_setup",*) 'Add timelabel?  : ', ATMOS_PHY_MP_RESTART_OUT_POSTFIX_TIMELABEL
    else
       LOG_INFO("ATMOS_PHY_MP_vars_setup",*) 'Restart output? : NO'
       ATMOS_PHY_MP_RESTART_OUTPUT = .false.
    endif


    ! diagnostices
    allocate( ATMOS_PHY_MP_CLDFRAC(KA,IA,JA,ADM_lall) )
    allocate( ATMOS_PHY_MP_Re     (KA,IA,JA,N_HYD,ADM_lall) )
    allocate( ATMOS_PHY_MP_Qe     (KA,IA,JA,N_HYD,ADM_lall) )
!OCL XFILL
    ATMOS_PHY_MP_CLDFRAC(:,:,:,:) = UNDEF
!OCL XFILL
    ATMOS_PHY_MP_Re     (:,:,:,:,:) = UNDEF
!OCL XFILL
    ATMOS_PHY_MP_Qe     (:,:,:,:,:) = UNDEF
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
  !> Open restart file for read
  subroutine ATMOS_PHY_MP_vars_restart_open

    return
  end subroutine ATMOS_PHY_MP_vars_restart_open

  !-----------------------------------------------------------------------------
  !> Read restart
  subroutine ATMOS_PHY_MP_vars_restart_read

    return
  end subroutine ATMOS_PHY_MP_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Create restart file
  subroutine ATMOS_PHY_MP_vars_restart_create

    return
  end subroutine ATMOS_PHY_MP_vars_restart_create

  !-----------------------------------------------------------------------------
  !> Exit netCDF define mode
  subroutine ATMOS_PHY_MP_vars_restart_enddef

    return
  end subroutine ATMOS_PHY_MP_vars_restart_enddef

  !-----------------------------------------------------------------------------
  !> Close restart file
  subroutine ATMOS_PHY_MP_vars_restart_close

    return
  end subroutine ATMOS_PHY_MP_vars_restart_close

  !-----------------------------------------------------------------------------
  !> Define variables in restart file
  subroutine ATMOS_PHY_MP_vars_restart_def_var

    return
  end subroutine ATMOS_PHY_MP_vars_restart_def_var

  !-----------------------------------------------------------------------------
  !> Write restart
  subroutine ATMOS_PHY_MP_vars_restart_write

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

    real(RP), intent(in) :: DENS(KA,IA,JA,ADM_lall)
    real(RP), intent(in) :: TEMP(KA,IA,JA,ADM_lall)
    real(RP), intent(in) :: QTRC(KA,IA,JA,QA,ADM_lall)

    real(RP) :: WORK  (KA,IA,JA,ADM_lall,N_HYD)
    logical  :: do_put
    integer  :: ih
    !---------------------------------------------------------------------------

    if ( HIST_CLDFRAC_id > 0 ) then
       call FILE_HISTORY_query( HIST_CLDFRAC_id, do_put )

       if ( do_put ) then
          call ATMOS_PHY_MP_vars_get_diagnostic( &
               DENS(:,:,:,:), TEMP(:,:,:,:), QTRC(:,:,:,:,:), & ! [IN]
               CLDFRAC=WORK(:,:,:,:,1)                        ) ! [OUT]
          call FILE_HISTORY_put( HIST_CLDFRAC_id, WORK(:,:,:,:,1) )
       end if
    end if

    if ( HIST_Re ) then
       do ih = 1, N_HYD
          if ( HIST_Re_id(ih) > 0 ) then
             call FILE_HISTORY_query( HIST_Re_id(ih), do_put )
             if ( do_put ) then
                call ATMOS_PHY_MP_vars_get_diagnostic( &
                     DENS(:,:,:,:), TEMP(:,:,:,:), QTRC(:,:,:,:,:), & ! [IN]
                     Re=WORK(:,:,:,:,:)                            ) ! [OUT]
                exit
             end if
          end if
       end do
       if ( do_put ) then
          do ih = 1, N_HYD
             if ( HIST_Re_id(ih) > 0 ) then
                call FILE_HISTORY_query( HIST_Re_id(ih), do_put )
                if ( do_put ) call FILE_HISTORY_put( HIST_Re_id(ih), WORK(:,:,:,:,ih) )
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
       ATMOS_PHY_MP_KESSLER_qtrc2qhyd, &
       ATMOS_PHY_MP_KESSLER_effective_radius, &
       ATMOS_PHY_MP_KESSLER_cloud_fraction
    use scale_atmos_phy_mp_tomita08, only: &
       ATMOS_PHY_MP_TOMITA08_qtrc2qhyd, &
       ATMOS_PHY_MP_TOMITA08_effective_radius, &
       ATMOS_PHY_MP_TOMITA08_cloud_fraction
    use scale_atmos_phy_mp_sn14, only: &
       ATMOS_PHY_MP_SN14_qtrc2qhyd, &
       ATMOS_PHY_MP_SN14_effective_radius, &
       ATMOS_PHY_MP_SN14_cloud_fraction
    use scale_atmos_phy_mp_suzuki10, only: &
       ATMOS_PHY_MP_suzuki10_qtrc2qhyd, &
       ATMOS_PHY_MP_suzuki10_effective_radius, &
       ATMOS_PHY_MP_suzuki10_cloud_fraction
    use mod_atmos_admin, only: &
       ATMOS_PHY_MP_TYPE

    real(RP), intent(in) :: DENS(KA,IA,JA,ADM_lall)
    real(RP), intent(in) :: TEMP(KA,IA,JA,ADM_lall)
    real(RP), intent(in) :: QTRC(KA,IA,JA,QA,ADM_lall)
    real(RP), intent(out), optional :: CLDFRAC(KA,IA,JA,ADM_lall)       !> cloud fraction [0-1]
    real(RP), intent(out), optional :: Re     (KA,IA,JA,ADM_lall,N_HYD) !> effective radius [cm]
    real(RP), intent(out), optional :: Qe     (KA,IA,JA,ADM_lall,N_HYD) !> mass ratio [kg/kg]

    integer :: k, i, j, l, ih

    if ( present(CLDFRAC) ) then
       if ( .not. DIAG_CLDFRAC ) then
          select case ( ATMOS_PHY_MP_TYPE )
          case ( 'KESSLER' )
             do l = 1, ADM_lall
                call ATMOS_PHY_MP_kessler_cloud_fraction( &
                     KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                     QTRC(:,:,:,QHS:QHE,l), ATMOS_PHY_MP_cldfrac_thleshold, & ! [IN]
                     ATMOS_PHY_MP_CLDFRAC(:,:,:,l)                          ) ! [OUT]
             end do
          case ( 'TOMITA08' )
             do l = 1, ADM_lall
                call ATMOS_PHY_MP_tomita08_cloud_fraction( &
                     KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                     QTRC(:,:,:,QHS:QHE,l), ATMOS_PHY_MP_cldfrac_thleshold, & ! [IN]
                     ATMOS_PHY_MP_CLDFRAC(:,:,:,l)                          ) ! [OUT]
             end do
          case ( 'SN14' )
             do l = 1, ADM_lall
                call ATMOS_PHY_MP_sn14_cloud_fraction( &
                     KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                     QTRC(:,:,:,QHS:QHE,l), ATMOS_PHY_MP_cldfrac_thleshold, & ! [IN]
                     ATMOS_PHY_MP_CLDFRAC(:,:,:,l)                          ) ! [OUT]
             end do
          case ( 'SUZUKI10' )
             do l = 1, ADM_lall
                call ATMOS_PHY_MP_suzuki10_cloud_fraction( &
                     KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                     QTRC(:,:,:,QHS:QHE,l), ATMOS_PHY_MP_cldfrac_thleshold, & ! [IN]
                     ATMOS_PHY_MP_CLDFRAC(:,:,:,l)                          ) ! [OUT]
             end do
          case default
!OCL XFILL
             ATMOS_PHY_MP_CLDFRAC(:,:,:,:) = 0.0_RP
          end select
          DIAG_CLDFRAC = .true.
       end if
!OCL XFILL
       do l = 1, ADM_lall
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          CLDFRAC(k,i,j,l) = ATMOS_PHY_MP_CLDFRAC(k,i,j,l)
       end do
       end do
       end do
       end do
    end if

    if ( present(Re) ) then
       if ( .not. DIAG_Re ) then
          select case ( ATMOS_PHY_MP_TYPE )
          case ( 'KESSLER' )
             do l = 1, ADM_lall
                call ATMOS_PHY_MP_kessler_effective_radius( &
                     KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                     DENS(:,:,:,l), TEMP(:,:,:,l), QTRC(:,:,:,QHS:QHE,l), & ! [IN]
                     ATMOS_PHY_MP_Re(:,:,:,:,l)                           ) ! [OUT]
             end do
          case ( 'TOMITA08' )
             do l = 1, ADM_lall
                call ATMOS_PHY_MP_tomita08_effective_radius( &
                     KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                     DENS(:,:,:,l), TEMP(:,:,:,l), QTRC(:,:,:,QHS:QHE,l), & ! [IN]
                     ATMOS_PHY_MP_Re(:,:,:,:,l)                           ) ! [OUT]
             end do
          case ( 'SN14' )
             do l = 1, ADM_lall
                call ATMOS_PHY_MP_sn14_effective_radius( &
                     KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                     DENS(:,:,:,l), TEMP(:,:,:,l), QTRC(:,:,:,QHS:QHE,l), & ! [IN]
                     ATMOS_PHY_MP_Re(:,:,:,:,l)                           ) ! [OUT]
             end do
          case ( 'SUZUKI10' )
             do l = 1, ADM_lall
                call ATMOS_PHY_MP_suzuki10_effective_radius( &
                     KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                     DENS(:,:,:,l), TEMP(:,:,:,l), QTRC(:,:,:,QHS:QHE,l), & ! [IN]
                     ATMOS_PHY_MP_Re(:,:,:,:,l)                           ) ! [OUT]
             end do
          case default
!OCL XFILL
             ATMOS_PHY_MP_Re(:,:,:,:,:) = 0.0_RP
          end select
          DIAG_Re = .true.
       end if
!OCL XFILL
       do ih = 1, N_HYD
       do l = 1, ADM_lall
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          Re(k,i,j,l,ih) = ATMOS_PHY_MP_Re(k,i,j,ih,l)
       end do
       end do
       end do
       end do
       end do
    end if

    if ( present(Qe) ) then
       if ( .not. DIAG_Qe ) then
          select case ( ATMOS_PHY_MP_TYPE )
          case ( 'KESSLER' )
             do l = 1, ADM_lall
                call ATMOS_PHY_MP_kessler_qtrc2qhyd( &
                     KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                     QTRC(:,:,:,QHS:QHE,l),     & ! [IN]
                     ATMOS_PHY_MP_Qe(:,:,:,:,l) ) ! [OUT]
             end do
          case ( 'TOMITA08' )
             do l = 1, ADM_lall
                call ATMOS_PHY_MP_tomita08_qtrc2qhyd( &
                     KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                     QTRC(:,:,:,QHS:QHE,l),     & ! [IN]
                     ATMOS_PHY_MP_Qe(:,:,:,:,l) ) ! [OUT]
             end do
          case ( 'SN14' )
             do l = 1, ADM_lall
                call ATMOS_PHY_MP_sn14_qtrc2qhyd( &
                     KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                     QTRC(:,:,:,QHS:QHE,l),     & ! [IN]
                     ATMOS_PHY_MP_Qe(:,:,:,:,l) ) ! [OUT]
             end do
          case ( 'SUZUKI10' )
             do l = 1, ADM_lall
                call ATMOS_PHY_MP_suzuki10_qtrc2qhyd( &
                     KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                     QTRC(:,:,:,QHS:QHE,l),     & ! [IN]
                     ATMOS_PHY_MP_Qe(:,:,:,:,l) ) ! [OUT]
             end do
          case default
!OCL XIFLL
             ATMOS_PHY_MP_Qe(:,:,:,:,:) = 0.0_RP
          end select
          DIAG_Qe = .true.
       end if
!OCL XIFLL
       do ih = 1, N_HYD
       do l = 1, ADM_lall
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          Qe(k,i,j,l,ih) = ATMOS_PHY_MP_Qe(k,i,j,ih,l)
       end do
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
