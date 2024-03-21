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
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: USER_tracer_setup
  public :: USER_setup
  public :: USER_finalize
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
  integer,                private, parameter :: QA = 1
  character(len=H_SHORT), private            :: QNAME(QA)
  character(len=H_MID),   private            :: QDESC(QA)
  character(len=H_SHORT), private            :: QUNIT(QA)

  data QNAME / 'QV' /

  data QDESC / 'Ratio of Water Vapor mass to total mass (Specific humidity)' /

  data QUNIT / 'kg/kg' /

  real(DP) :: t_npf = 21600.D0 ! duration time for new particle formation

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Tracer setup
  subroutine USER_tracer_setup
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_regist
    use mod_atmos_phy_mp_vars, only: &
       QA_MP, &
       QS_MP, &
       QE_MP
    use mod_atmos_phy_mp_driver, only: &
       ATMOS_PHY_MP_USER_qhyd2qtrc
    implicit none

    !---------------------------------------------------------------------------

    call ATMOS_HYDROMETEOR_regist( 0, 0,                & ! (in)
                                   QNAME, QDESC, QUNIT, & ! (in)
                                   QS_MP                ) ! (out)
    QA_MP = 1
    QE_MP = QS_MP

    ATMOS_PHY_MP_USER_qhyd2qtrc => USER_qhyd2qtrc

    return
  end subroutine USER_tracer_setup

  !-----------------------------------------------------------------------------
  !> Setup before setup of other components
  subroutine USER_setup
    use scale_prc, only: &
       PRC_abort
    implicit none

    namelist / PARAM_USER / &
       t_npf

    integer :: ierr
    !---------------------------------------------------------------------------

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
  !> Finalization
  subroutine USER_finalize
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_finalize

  !-----------------------------------------------------------------------------
  !> Make initial state
  subroutine USER_mkinit
    implicit none
    !---------------------------------------------------------------------------

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
  !> User step
  subroutine USER_update
    use scale_time, only: &
       TIME_NOWDAYSEC, &
       TIME_STARTDAYSEC
    use scale_atmos_phy_ae_kajino13, only: &
       ATMOS_PHY_AE_KAJINO13_flag_npf
    implicit none
    !---------------------------------------------------------------------------
    real(DP) :: t_elaps

    t_elaps = TIME_NOWDAYSEC - TIME_STARTDAYSEC

    if ( t_elaps > t_npf ) then ! no more new particle formation does not occur
       ATMOS_PHY_AE_KAJINO13_flag_npf = .false.
    end if

    return
  end subroutine USER_update

  subroutine USER_qhyd2qtrc( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       QV, QHYD, &
       QTRC, &
       QNUM  )
    use scale_atmos_hydrometeor, only: &
         N_HYD
    use mod_atmos_phy_mp_vars, only: &
         QA_MP
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in) :: QV   (KA,IA,JA)
    real(RP), intent(in) :: QHYD(KA,IA,JA,N_HYD)

    real(RP), intent(out) :: QTRC(KA,IA,JA,QA_MP)

    real(RP), intent(in), optional :: QNUM(KA,IA,JA,N_HYD)

    QTRC(:,:,:,1) = QV(:,:,:)

    return
  end subroutine USER_qhyd2qtrc

end module mod_user
