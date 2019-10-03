!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User defined module
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_user
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
  public :: USER_tracer_setup
  public :: USER_setup
  public :: USER_mkinit
  public :: USER_calc_tendency
  public :: USER_update

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
  logical, private :: USER_do = .false. !< do user step?
  integer, private, parameter :: I_num = 2
  real(RP), private :: position_phi_x(I_num)  ![m]
  real(RP), private :: position_phi_y(I_num)  ![m]
  real(RP), private :: position_phi_z(I_num)  ![m]
  real(RP), private :: radius_phi_x(I_num) ![m]
  real(RP), private :: radius_phi_y(I_num) ![m]
  real(RP), private :: radius_phi_z(I_num) ![m]
  real(RP), private :: qvalue(I_num)   ![m]
  !--- Indeces for determining species of cloud particle
!  integer, parameter, private :: LTS = 0
  !
  integer,  private, parameter   :: HYDRO_MAX = 5
  integer,  private, parameter   :: I_QV = 1
  integer,  private, parameter   :: I_QC = 2
  integer,  private, parameter   :: I_QR = 3
  integer,  private, parameter   :: I_QI = 4
  integer,  private, parameter   :: I_QS = 5
  integer,  private, parameter   :: I_QG = 6
  integer,  private, parameter   :: I_NC = 7
  integer,  private, parameter   :: I_NR = 8
  integer,  private, parameter   :: I_NI = 9
  integer,  private, parameter   :: I_NS = 10
  integer,  private, parameter   :: I_NG = 11

  integer, private, parameter :: QA_MP_TOMITA08 = 6
  integer, private, parameter :: QA_MP_SN14 = 11
  integer, private, parameter :: I_mp_QC = 1
  integer, private, parameter :: I_mp_QR = 2
  integer, private, parameter :: I_mp_QI = 3
  integer, private, parameter :: I_mp_QS = 4
  integer, private, parameter :: I_mp_QG = 5
  integer, private, parameter :: I_mp_NC = 6
  integer, private, parameter :: I_mp_NR = 7
  integer, private, parameter :: I_mp_NI = 8
  integer, private, parameter :: I_mp_NS = 9
  integer, private, parameter :: I_mp_NG = 10
  real(RP), private :: q_rate(HYDRO_MAX) = (/1.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP /)
  real(RP), private :: q_typ(I_mp_QC:I_mp_QG) = (/0.1_RP, 0.1_RP, 0.1_RP, 0.1_RP, 0.1_RP /) ! [g/kg]
  real(RP), private :: n_typ(I_mp_NC:I_mp_NG) = (/100.0_RP, 1.0_RP, 10.0_RP, 10.0_RP, 0.1_RP /) ! [/cm3]
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Config before setup of tracers
  subroutine USER_tracer_setup
    use scale_tracer, only: &
       TRACER_regist
    implicit none

    ! if you want to add tracers, call the TRACER_regist subroutine.
    ! e.g.,
!    integer, parameter     :: NQ = 1
!    integer                :: QS
!    character(len=H_SHORT) :: NAME(NQ)
!    character(len=H_MID)   :: DESC(NQ)
!    character(len=H_SHORT) :: UNIT(NQ)
!
!    data NAME (/ 'name' /)
!    data DESC (/ 'tracer name' /)
!    data UNIT (/ 'kg/kg' /)
    !---------------------------------------------------------------------------

!    call TRACER_regist( QS,   & ! [OUT]
!                        NQ,   & ! [IN]
!                        NAME, & ! [IN]
!                        DESC, & ! [IN]
!                        UNIT  ) ! [IN]

    return
  end subroutine USER_tracer_setup

  !-----------------------------------------------------------------------------
  !> Setup before setup of other components
  subroutine USER_setup
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       PI    => CONST_PI
    use scale_atmos_grid_cartesC, only: &
       CDX   => ATMOS_GRID_CARTESC_CDX, &
       CDY   => ATMOS_GRID_CARTESC_CDY, &
       CX    => ATMOS_GRID_CARTESC_CX, &
       CY    => ATMOS_GRID_CARTESC_CY, &
       CZ    => ATMOS_GRID_CARTESC_CZ
    use scale_file_history, only: &
       FILE_HISTORY_reg
    use mod_atmos_admin, only: &
       ATMOS_sw_phy_mp, &
       ATMOS_PHY_MP_TYPE
    implicit none

    namelist / PARAM_USER / &
       USER_do, &
       position_phi_x, &
       position_phi_y, &
       position_phi_z, &
       radius_phi_x, &
       radius_phi_y, &
       radius_phi_z, &
       qvalue, &
       q_rate, &
       q_typ, &
       n_typ

    integer :: ierr, ip
    !---------------------------------------------------------------------------

    if( ATMOS_PHY_MP_TYPE /= 'TOMITA08' .and. ATMOS_PHY_MP_TYPE /= 'SN14' ) then
       LOG_ERROR("USER_setup",*) 'MP_TYPE should be TOMITA08 or SN14 for this test. Check!'
       call PRC_abort
    endif

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_USER,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("USER_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("USER_setup",*) 'Not appropriate names in namelist PARAM_USER. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_USER)

    LOG_NEWLINE
    LOG_INFO("USER_setup",*) 'Setup'

    do ip = 1, I_num
       if(  radius_phi_x(ip) /= radius_phi_y(ip) .or. &
            radius_phi_y(ip) /= radius_phi_z(ip) .or. &
            radius_phi_x(ip) /= radius_phi_z(ip)      ) then
          LOG_ERROR("USER_setup",*) 'xxx radius_phi_x, _y, and _z should be same, stop!', ip
          call PRC_abort
       endif
    enddo

    ATMOS_sw_phy_mp = .false.

    return
  end subroutine USER_setup

  !-----------------------------------------------------------------------------
  !> Make initial state
  subroutine USER_mkinit
    use mod_atmos_vars, only: &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       RHOT, &
       QTRC, &
       POTT
    use scale_atmos_grid_cartesC, only: &
       CX    => ATMOS_GRID_CARTESC_CX, &
       CY    => ATMOS_GRID_CARTESC_CY, &
       CZ    => ATMOS_GRID_CARTESC_CZ
    use scale_const, only: &
       UNDEF => CONST_UNDEF, &
       SMALL => CONST_EPS
    use mod_atmos_admin, only: &
       ATMOS_PHY_MP_TYPE
    use mod_atmos_phy_lt_vars, only: &
       ATMOS_PHY_LT_Epot, &
       QS_LT, QE_LT
    implicit none
    real(RP) :: distance(3), total_rate, qc, qr, qi, qs, qg, nc, nr, ni, ns, ng
    integer  :: k, i, j, ip, iq
    !---------------------------------------------------------------------------

    LOG_INFO("USER_mkinit",*) 'Setup USER_mkinit'
    QTRC(:,:,:,:) = 0.0_RP
    DENS(:,:,:) = 1.0_RP
    MOMX(:,:,:) = 0.0_RP
    MOMY(:,:,:) = 0.0_RP
    MOMZ(:,:,:) = 0.0_RP
    pott(:,:,:) = 300.0_RP
    RHOT(:,:,:) = DENS(:,:,:) * pott(:,:,:)

    NC = UNDEF
    select case ( ATMOS_PHY_MP_TYPE )
    case ( 'TOMITA08' )
       qc = q_typ(I_mp_QC)*1.0E-3_RP  ![g/kg] -> [kg/kg]
       qr = q_typ(I_mp_QR)*1.0E-3_RP  ![g/kg] -> [kg/kg]
       qi = q_typ(I_mp_QI)*1.0E-3_RP  ![g/kg] -> [kg/kg]
       qs = q_typ(I_mp_QS)*1.0E-3_RP  ![g/kg] -> [kg/kg]
       qg = q_typ(I_mp_QG)*1.0E-3_RP  ![g/kg] -> [kg/kg]
    case ( 'SN14' )
       qc = q_typ(I_mp_QC)*1.0E-3_RP  ![g/kg] -> [kg/kg]
       qr = q_typ(I_mp_QR)*1.0E-3_RP  ![g/kg] -> [kg/kg]
       qi = q_typ(I_mp_QI)*1.0E-3_RP  ![g/kg] -> [kg/kg]
       qs = q_typ(I_mp_QS)*1.0E-3_RP  ![g/kg] -> [kg/kg]
       qg = q_typ(I_mp_QG)*1.0E-3_RP  ![g/kg] -> [kg/kg]
       NC = n_typ(I_mp_NC)*1.0E+6_RP  ! [/cm3] -> [/m3]
       NR = n_typ(I_mp_NR)*1.0E+6_RP  ! [/cm3] -> [/m3]
       NI = n_typ(I_mp_NI)*1.0E+6_RP  ! [/cm3] -> [/m3]
       NS = n_typ(I_mp_NS)*1.0E+6_RP  ! [/cm3] -> [/m3]
       NG = n_typ(I_mp_NG)*1.0E+6_RP  ! [/cm3] -> [/m3]
    end select

    total_rate = 0.0_RP
    do iq = 1, HYDRO_MAX
       total_rate = total_rate + q_rate(iq)*( 0.5_RP + sign( 0.50_RP, q_typ(iq+1)-SMALL ) )
    enddo

    do ip = 1, I_num
    do j = 1, JA-1
    do i = 1, IA-1
    do k = 1, KA-1
       distance(3) = sqrt( ( CZ(k)-position_phi_z(ip) )**2 &
                         + ( CX(i)-position_phi_x(ip) )**2 &
                         + ( CY(j)-position_phi_y(ip) )**2 )

       if( distance(3) < radius_phi_z(ip) ) then
          do iq = QS_LT, QE_LT
            QTRC(k,i,j,iq) = qvalue(ip) * 1.0E+15_RP &
                           * q_rate(iq-QS_LT+1)/total_rate &
                           * ( 0.5_RP + sign( 0.50_RP, q_typ(iq-QS_LT+1)-SMALL ) )![fC/kg]
          enddo
          QTRC(k,i,j,I_QC) = qc
          QTRC(k,i,j,I_QR) = qr
          QTRC(k,i,j,I_QI) = qi
          QTRC(k,i,j,I_QS) = qs
          QTRC(k,i,j,I_QG) = qg
          if( NC /= UNDEF ) then
             QTRC(k,i,j,I_NC) = nc
             QTRC(k,i,j,I_NR) = nr
             QTRC(k,i,j,I_NI) = ni
             QTRC(k,i,j,I_NS) = ns
             QTRC(k,i,j,I_NG) = ng
          endif
       endif

    enddo
    enddo
    enddo
    enddo

    ATMOS_PHY_LT_Epot (:,:,:)   = 0.0_RP

    return
  end subroutine USER_mkinit

  !-----------------------------------------------------------------------------
  !> Calculation tendency
  subroutine USER_calc_tendency
    implicit none

    return
  end subroutine USER_calc_tendency

  !-----------------------------------------------------------------------------
  !> User step
  !> User step
  subroutine USER_update
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_phy_mp_tomita08, only: &
       ATMOS_PHY_MP_tomita08_Sarea
    use scale_atmos_phy_mp_sn14, only: &
       ATMOS_PHY_MP_sn14_Sarea
    use scale_atmos_phy_mp_suzuki10, only: &
       ATMOS_PHY_MP_suzuki10_Sarea
    use mod_atmos_phy_lt_vars, only: &
       Sarea => ATMOS_PHY_LT_Sarea
    use mod_atmos_admin, only: &
       ATMOS_PHY_MP_TYPE
    implicit none
    integer :: i, j, k
    real(RP) :: dummy_dens(KA,IA,JA)
    real(RP), allocatable :: dummy_qtrc(:,:)
    !---------------------------------------------------------------------------

    if ( USER_do ) then
       dummy_dens(:,:,:) = 1.0_RP 
       select case ( ATMOS_PHY_MP_TYPE )
       case ( 'TOMITA08' )

          allocate(dummy_qtrc(KA,QA_MP_TOMITA08))
          dummy_qtrc(:,I_QC) = q_typ(I_mp_QC)*1.0E-3_RP  ![g/kg] -> [kg/kg]
          dummy_qtrc(:,I_QR) = q_typ(I_mp_QR)*1.0E-3_RP  ![g/kg] -> [kg/kg]
          dummy_qtrc(:,I_QI) = q_typ(I_mp_QI)*1.0E-3_RP  ![g/kg] -> [kg/kg]
          dummy_qtrc(:,I_QS) = q_typ(I_mp_QS)*1.0E-3_RP  ![g/kg] -> [kg/kg]
          dummy_qtrc(:,I_QG) = q_typ(I_mp_QG)*1.0E-3_RP  ![g/kg] -> [kg/kg]
          do j = 1, JA
          do i = 1, IA
             call ATMOS_PHY_MP_tomita08_Sarea( KA, KS, KE, &
                                               dummy_dens(:,i,j), &
                                               dummy_qtrc(:,:),   &
                                               Sarea(:,i,j,:)     )
          enddo
          enddo
          deallocate(dummy_qtrc)

       case ( 'SN14' )

          allocate(dummy_qtrc(KA,QA_MP_SN14))
          dummy_qtrc(:,I_QC) = q_typ(I_mp_QC)*1.0E-3_RP  ! [g/kg] -> [kg/kg]
          dummy_qtrc(:,I_QR) = q_typ(I_mp_QR)*1.0E-3_RP  ! [g/kg] -> [kg/kg]
          dummy_qtrc(:,I_QI) = q_typ(I_mp_QI)*1.0E-3_RP  ! [g/kg] -> [kg/kg]
          dummy_qtrc(:,I_QS) = q_typ(I_mp_QS)*1.0E-3_RP  ! [g/kg] -> [kg/kg]
          dummy_qtrc(:,I_QG) = q_typ(I_mp_QG)*1.0E-3_RP  ! [g/kg] -> [kg/kg]
          dummy_qtrc(:,I_NC) = n_typ(I_mp_NC)*1.0E+6_RP  ! [/cm3] -> [/m3]
          dummy_qtrc(:,I_NR) = n_typ(I_mp_NR)*1.0E+6_RP  ! [/cm3] -> [/m3]
          dummy_qtrc(:,I_NI) = n_typ(I_mp_NI)*1.0E+6_RP  ! [/cm3] -> [/m3]
          dummy_qtrc(:,I_NS) = n_typ(I_mp_NS)*1.0E+6_RP  ! [/cm3] -> [/m3]
          dummy_qtrc(:,I_NG) = n_typ(I_mp_NG)*1.0E+6_RP  ! [/cm3] -> [/m3]
          do j = 1, JA
          do i = 1, IA
             call ATMOS_PHY_MP_sn14_Sarea( KA, KS, KE, &
                                           dummy_dens(:,i,j), &
                                           dummy_qtrc(:,:),   &
                                           Sarea(:,i,j,:)     )
          enddo
          enddo
          deallocate(dummy_qtrc)

!       case ( 'SUZUKI10' )
!          do j = 1, JA
!          do i = 1, IA
!             call ATMOS_PHY_MP_suzuki10_Sarea( KA, KS, KE, &
!                                               dummy_dens(:,i,j), &
!                                               dummy_qtrc(:,:),   &
!                                               Sarea(:,i,j,:)     )
!          enddo
!          enddo
       end select

    endif
    LOG_INFO("USER_update",*) 'Sarea', Sarea(KS,IS,JS,:)

    return
  end subroutine USER_update

end module mod_user
