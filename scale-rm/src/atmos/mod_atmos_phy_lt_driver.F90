!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Chemistry
!!
!! @par Description
!!          Atmospheric lightning driver
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_atmos_phy_lt_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
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
  public :: ATMOS_PHY_LT_driver_tracer_setup
  public :: ATMOS_PHY_LT_driver_setup
  public :: ATMOS_PHY_LT_driver_finalize
  public :: ATMOS_PHY_LT_driver_adjustment

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !

  !--- For history output
  integer, private, parameter :: w_nmax = 3
  integer, private, parameter :: I_CRGD_LIQ = 1
  integer, private, parameter :: I_CRGD_ICE = 2
  integer, private, parameter :: I_CRGD_TOT = 3
  integer,  private              :: HIST_id(w_nmax)
  character(len=H_SHORT), private :: w_name(w_nmax)
  character(len=H_MID),   private :: w_longname(w_nmax)
  character(len=H_SHORT), private :: w_unit(w_nmax)
  data w_name / 'CRGD_LIQ', &
                'CRGD_ICE', &
                'CRGD_TOT' /
  data w_longname / 'Charge density of liquid water', &
                    'Charge density of ice water', &
                    'Charge density of QHYD' /
  data w_unit / 'nC/m3', &
                'nC/m3', &
                'nC/m3' /


  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Config
  subroutine ATMOS_PHY_LT_driver_tracer_setup
    use mod_atmos_admin, only: &
       ATMOS_PHY_MP_TYPE, &
       ATMOS_PHY_LT_TYPE, &
       ATMOS_sw_phy_lt
    use scale_tracer, only: &
       TRACER_regist
    use scale_atmos_phy_mp_suzuki10, only: &
       ATMOS_PHY_MP_suzuki10_nwaters, &
       ATMOS_PHY_MP_suzuki10_nices, &
       ATMOS_PHY_MP_suzuki10_nccn
    use mod_atmos_phy_lt_vars, only: &
       QA_LT, &
       QS_LT, &
       QE_LT
    use scale_atmos_hydrometeor, only: &
       QHA, &
       QHS
    use scale_prc, only: &
       PRC_abort
    implicit none

    character(len=H_SHORT), allocatable :: NAME(:)
    character(len=H_MID  ), allocatable :: DESC(:)
    character(len=H_SHORT), allocatable :: UNIT(:)
    integer :: iq
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_LT_driver_tracer_setup",*) 'Setup'

    if ( ATMOS_sw_phy_lt ) then

       select case ( ATMOS_PHY_LT_TYPE )
       case ( 'OFF', 'NONE' )
          LOG_INFO("ATMOS_PHY_LT_driver_tracer_setup",*) 'this component is never called.'
       case ( 'SATO2019' )
          select case ( ATMOS_PHY_MP_TYPE )
          case ( 'TOMITA08', 'SN14' )
             ! do nothing
          case ( 'SUZUKI10' )

             if( ATMOS_PHY_MP_suzuki10_nccn /= 0 ) then
                LOG_ERROR("ATMOS_PHY_LT_driver_tracer_setup",*) 'nccn in SUZUKI10 should be 0 for lithgning component(', ATMOS_PHY_MP_suzuki10_nccn, '). CHECK!'
                call PRC_abort
             endif
             if ( ATMOS_PHY_MP_suzuki10_nices == 0 ) then
                LOG_ERROR("ATMOS_PHY_LT_driver_tracer_setup",*) 'ICEFLG in SUZUKI10 should be 1 for lithgning component. CHECK!'
                call PRC_abort
             endif
          case ( 'KESSLER' )
             LOG_ERROR("ATMOS_PHY_LT_driver_tracer_setup",*) 'ATMOS_PHY_MP_TYPE should be TOMITA08, or SN14, or SUZUKI10 (', ATMOS_PHY_MP_TYPE, '). CHECK!'
             call PRC_abort
          end select
       case default
          LOG_ERROR("ATMOS_PHY_LT_driver_tracer_setup",*) 'invalid lithgning type(', ATMOS_PHY_LT_TYPE, '). CHECK!'
          call PRC_abort
       end select

       QA_LT = QHA

       allocate( NAME(QA_LT), DESC(QA_LT), UNIT(QA_LT) )
       do iq = 1, QA_LT
          NAME(iq) = 'QCRG_'//trim(TRACER_NAME(QHS+iq-1)(2:))
          DESC(iq) = 'Ratio of charge density of '//trim(TRACER_NAME(QHS+iq-1))
          UNIT(iq) = 'fC/kg'
       end do
       call TRACER_regist( QS_LT,   & ! [OUT]
                           QA_LT,   & ! [IN]
                           NAME(:), & ! [IN]
                           DESC(:), & ! [IN]
                           UNIT(:)  ) ! [IN]
       deallocate( NAME, DESC, UNIT )

       QE_LT = QS_LT - 1 + QA_LT

    else
       QA_LT = 0
       QS_LT = -1
       QE_LT = -2
    endif

    return
  end subroutine ATMOS_PHY_LT_driver_tracer_setup

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_LT_driver_setup
    use scale_atmos_grid_cartesC_real, only: &
       REAL_LON => ATMOS_GRID_CARTESC_REAL_LON, &
       REAL_LAT => ATMOS_GRID_CARTESC_REAL_LAT
    use scale_atmos_grid_cartesC, only: &
       CDX => ATMOS_GRID_CARTESC_CDX, &
       CDY => ATMOS_GRID_CARTESC_CDY
    use scale_atmos_phy_lt_sato2019, only: &
       ATMOS_PHY_LT_sato2019_setup
    use mod_atmos_phy_lt_vars, only: &
       flg_lt, &
       ATMOS_PHY_LT_Sarea
    use mod_atmos_admin, only: &
       ATMOS_PHY_MP_TYPE, &
       ATMOS_PHY_LT_TYPE
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_hydrometeor, only: &
       QHA
    use scale_file_history, only: &
       FILE_HISTORY_reg
    implicit none

    logical :: LT_force_with_SUZUKI10 = .false. ! experimental use only

    namelist / PARAM_ATMOS_PHY_LT / &
         LT_force_with_SUZUKI10

    integer :: ip
    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_LT_driver_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_LT,iostat=ierr)
    if ( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_PHY_LT_driver_setup",*) 'Not found namelist. Default used.'
    elseif ( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_LT_dirver_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_LT. Check!'
       call PRC_abort
    end if
    LOG_NML(PARAM_ATMOS_PHY_LT)

    select case( ATMOS_PHY_LT_TYPE )
    case ( 'SATO2019' )

       if ( ATMOS_PHY_MP_TYPE == "SUZUKI10" ) then
          LOG_WARN("ATMOS_PHY_LT_driver_setup",*) 'At this moment, ATMOS_PHY_MP_TYPE of SUZUKI10 with SATO2019 scheme is for experimental use only.'
          if ( .not. LT_force_with_SUZUKI10 ) then
             LOG_ERROR("ATMOS_PHY_LT_driver_setup",*) 'At this moment, ATMOS_PHY_MP_TYPE of SUZUKI10 with SATO2019 scheme is for experimental use only.'
             call PRC_abort
          end if
       end if


       call ATMOS_PHY_LT_sato2019_setup( KA, KS, KE, & ! [IN]
                                         IA, IS, IE, & ! [IN]
                                         JA, JS, JE, & ! [IN]
                                         IMAXG,      & ! [IN]
                                         JMAXG,      & ! [IN]
                                         KMAX,       & ! [IN]
                                         ATMOS_PHY_MP_TYPE, & ! [IN]
                                         CDX, CDY           ) ! [IN]
       flg_lt = .true.
    case default
       flg_lt = .false.
    end select


    if ( flg_lt ) then

       allocate( ATMOS_PHY_LT_Sarea(KA,IA,JA,QHA) )

       do ip = 1, w_nmax
          call FILE_HISTORY_reg( w_name(ip), w_longname(ip), w_unit(ip), & ! [IN]
                                 HIST_id(ip)                             ) ! [OUT]
       end do

       call history

    else
       LOG_INFO("ATMOS_PHY_LT_driver_setup",*) 'This component is never called.'
    endif


    return
  end subroutine ATMOS_PHY_LT_driver_setup

  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_LT_driver_finalize
    use scale_atmos_phy_lt_sato2019, only: &
       ATMOS_PHY_LT_sato2019_finalize
    use mod_atmos_admin, only: &
       ATMOS_PHY_LT_TYPE, &
       ATMOS_sw_phy_lt
    use mod_atmos_phy_lt_vars, only: &
       flg_lt, &
       ATMOS_PHY_LT_Sarea
    implicit none
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_LT_driver_finalize",*) 'Finalize'

    if ( ATMOS_sw_phy_lt ) then
       select case ( ATMOS_PHY_LT_TYPE )
       case ( 'OFF', 'NONE' )
       case ( 'SATO2019' )
          call ATMOS_PHY_LT_sato2019_finalize
       end select
    end if

    if ( flg_lt ) then
       deallocate( ATMOS_PHY_LT_Sarea )
    end if

    return
  end subroutine ATMOS_PHY_LT_driver_finalize
  !-----------------------------------------------------------------------------
  !> Driver
  subroutine ATMOS_PHY_LT_driver_adjustment
    use scale_time, only: &
       dt_LT => TIME_DTSEC_ATMOS_PHY_LT
    use mod_atmos_vars, only: &
       DENS => DENS_av, &
       RHOT => RHOT_av, &
       QTRC => QTRC_av, &
       ATMOS_vars_get_diagnostic
    use mod_atmos_phy_lt_vars, only: &
       QA_LT, &
       QS_LT, &
       QE_LT, &
       Sarea => ATMOS_PHY_LT_Sarea, &
       Epot  => ATMOS_PHY_LT_Epot
    use scale_atmos_phy_lt_sato2019, only: &
       ATMOS_PHY_LT_sato2019_adjustment
    implicit none

    real(RP) :: QHYD(KA,IA,JA)
    !---------------------------------------------------------------------------

    !$acc data create(QHYD) &
    !$acc      copyin(DENS,RHOT,QHYD,Sarea) &
    !$acc      copy(QTRC,Epot)

    call ATMOS_vars_get_diagnostic( "QHYD", QHYD(:,:,:) )

    call ATMOS_PHY_LT_sato2019_adjustment( &
         KA, KS, KE, IA, IS, IE, JA, JS, JE, KIJMAX, IMAX, JMAX, QA_LT, & ! [IN]
         DENS(:,:,:), RHOT(:,:,:), QHYD(:,:,:), Sarea(:,:,:,:), & ! [IN]
         dt_LT,                                                 & ! [IN]
         QTRC(:,:,:,QS_LT:QE_LT), Epot(:,:,:)                   ) ! [INOUT]

    call history

    !$acc end data

    return
  end subroutine ATMOS_PHY_LT_driver_adjustment

  ! private
  subroutine history
    use scale_atmos_hydrometeor, only: &
       QLA
    use mod_atmos_phy_lt_vars, only: &
       QS_LT, &
       QE_LT
    use mod_atmos_vars, only: &
       DENS => DENS_av, &
       QTRC => QTRC_av
    use scale_file_history, only: &
       FILE_HISTORY_query, &
       FILE_HISTORY_put
    implicit none
    real(RP) :: work(KA,IA,JA)
    logical  :: HIST_sw(w_nmax)
    integer  :: k, i, j, n, ip

    !$acc data copyin(QTRC,DENS) create(work)

    do ip = 1, w_nmax
       call FILE_HISTORY_query( HIST_id(ip), HIST_sw(ip) )
    end do

    if ( HIST_sw(I_CRGD_LIQ) ) then
       !$omp parallel do
       !$acc kernels
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          work(k,i,j) = 0.0_RP
          do n = QS_LT, QS_LT + QLA - 1
             work(k,i,j) = work(k,i,j) + QTRC(k,i,j,n)
          enddo
          work(k,i,j) = work(k,i,j) * DENS(k,i,j) * 1.E-6_RP ! [fC/kg] -> [nc/m3]
       end do
       end do
       end do
       !$acc end kernels
       call FILE_HISTORY_put( HIST_id(I_CRGD_LIQ), work(:,:,:) )
    end if
    if ( HIST_sw(I_CRGD_ICE) ) then
       !$omp parallel do
       !$acc kernels
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          work(k,i,j) = 0.0_RP
          do n = QS_LT + QLA, QE_LT
             work(k,i,j) = work(k,i,j) + QTRC(k,i,j,n)
          enddo
          work(k,i,j) = work(k,i,j) * DENS(k,i,j) * 1.E-6_RP ! [fC/kg] -> [nc/m3]
       end do
       end do
       end do
       !$acc end kernels
       call FILE_HISTORY_put( HIST_id(I_CRGD_ICE), work(:,:,:) )
    end if
    if ( HIST_sw(I_CRGD_TOT) ) then
       !$omp parallel do
       !$acc kernels
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          work(k,i,j) = 0.0_RP
          do n = QS_LT, QE_LT
             work(k,i,j) = work(k,i,j) + QTRC(k,i,j,n)
          enddo
          work(k,i,j) = work(k,i,j) * DENS(k,i,j) * 1.E-6_RP ! [fC/kg] -> [nc/m3]
       end do
       end do
       end do
       !$acc end kernels
       call FILE_HISTORY_put( HIST_id(I_CRGD_TOT), work(:,:,:) )
    end if

    !$acc end data

    return
  end subroutine history

end module mod_atmos_phy_lt_driver
