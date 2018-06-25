!-------------------------------------------------------------------------------
!> module User
!!
!! @par Description
!!          calc perturbation
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
  logical,  private, save :: USER_do = .false. !< do user step?

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Tracer setup
  subroutine USER_tracer_setup

    return
  end subroutine USER_tracer_setup

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine USER_setup
    use scale_prc, only: &
       PRC_abort
    implicit none

    namelist / PARAM_USER / &
       USER_do

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("USER_setup",*) 'Setup'
    LOG_INFO("USER_setup",*) 'User procedure in test/framework/terrainfollowing/500m_precip'

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

    return
  end subroutine USER_setup

  !-----------------------------------------------------------------------------
  !> Make initial state
  subroutine USER_mkinit
    use scale_atmos_grid_cartesC_real, only : &
       CZ => ATMOS_GRID_CARTESC_REAL_CZ, &
       FZ => ATMOS_GRID_CARTESC_REAL_FZ
    use scale_atmos_hydrometeor, only: &
       N_HYD, &
       I_HR
    use mod_atmos_vars, only: &
       QTRC
    use mod_atmos_phy_mp_driver, only: &
       ATMOS_PHY_MP_driver_qhyd2qtrc
    use mod_atmos_phy_mp_vars, only: &
       QS_MP, &
       QE_MP
    implicit none

    real(RP) :: QV  (KA,IA,JA)
    real(RP) :: QHYD(KA,IA,JA,N_HYD)
    real(RP) :: QNUM(KA,IA,JA,N_HYD)

    integer  :: modsec
    real(RP) :: dist

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    LOG_INFO("USER_mkinit",*) 'Add rain.'

    QV  (:,:,:)   = 0.0_RP
    QHYD(:,:,:,:) = 0.0_RP
    QNUM(:,:,:,:) = 0.0_RP

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
!        dist = ( ( CZ(k,i,j) - 3000.0_RP ) / 50.0_RP )**2
!
!        QHYD(k,i,j,I_HR) = 1.E-3 / ( 1.0_RP + dist )
       if (       FZ(k,  i,j) >= 3000.0_RP &
            .AND. FZ(k-1,i,j) <  3000.0_RP ) then
          QHYD(k,i,j,I_HR) = 1.E-3_RP
       endif
    enddo
    enddo
    enddo

    call ATMOS_PHY_MP_driver_qhyd2qtrc( KA, KS, KE,              & ! [IN]
                                        IA, IS, IE,              & ! [IN]
                                        JA, JS, JE,              & ! [IN]
                                        QV  (:,:,:),             & ! [IN]
                                        QHYD(:,:,:,:),           & ! [IN]
                                        QTRC(:,:,:,QS_MP:QE_MP), & ! [OUT]
                                        QNUM=QNUM(:,:,:,:)       ) ! [IN]

    return
  end subroutine USER_mkinit

  !-----------------------------------------------------------------------------
  !> Calculate tendency
  subroutine USER_calc_tendency
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_calc_tendency

  !-----------------------------------------------------------------------------
  !> Step
  subroutine USER_update
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_update

end module mod_user
