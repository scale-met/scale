!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User defined module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2012-12-26 (H.Yashiro)   [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_user
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
  public :: USER_config
  public :: USER_setup
  public :: USER_resume0
  public :: USER_resume
  public :: USER_step

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
  integer :: I_NC

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Config
  subroutine USER_config
    use scale_tracer, only: &
       TRACER_regist
    use mod_atmos_phy_mp_vars, only: &
         QA_MP, &
         QS_MP, &
         QE_MP
    use mod_atmos_phy_mp_driver, only: &
         ATMOS_PHY_MP_USER_qhyd2qtrc
    implicit none
    !---------------------------------------------------------------------------

    call TRACER_REGIST( QS_MP,                & ! [OUT]
                        1,                    & ! [IN]
                        (/'NC'/),             & ! [IN]
                        (/'Passive tracer'/), & ! [IN]
                        (/'1'/)               ) ! [IN]

    QA_MP = 1
    QE_MP = QS_MP

    ATMOS_PHY_MP_USER_qhyd2qtrc => USER_qhyd2qtrc

    return
  end subroutine USER_config

  !-----------------------------------------------------------------------------
  !> Setup before setup of other components
  subroutine USER_setup
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_setup

  !-----------------------------------------------------------------------------
  !> Resuming operation, before calculating tendency
  subroutine USER_resume0
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_resume0

  !-----------------------------------------------------------------------------
  !> Resuming operation
  subroutine USER_resume
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_resume

  !-----------------------------------------------------------------------------
  !> User step
  subroutine USER_step
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_step

  subroutine USER_qhyd2qtrc( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       QV, QHYD, &
       QTRC, &
       QNUM  )
    use scale_atmos_hydrometeor, only: &
         N_HYD, &
         I_HC
    use mod_atmos_phy_mp_vars, only: &
         QA_MP
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in) :: QV   (KA,IA,JA)
    real(RP), intent(in) :: QHYD(KA,IA,JA,N_HYD)

    real(RP), intent(out) :: QTRC(KA,IA,JA,QA_MP)

    real(RP), intent(in), optional :: QNUM(KA,IA,JA,N_HYD)

    QTRC(:,:,:,1) = QNUM(:,:,:,I_HC)

    return
  end subroutine USER_qhyd2qtrc

end module mod_user
