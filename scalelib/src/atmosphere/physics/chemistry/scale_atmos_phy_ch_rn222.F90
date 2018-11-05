!-------------------------------------------------------------------------------
!> module atmosphere / physics / chemistry / RN222
!!
!! @par Description
!!         Component for rn222 tracer
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_phy_ch_rn222
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_CH_rn222_setup
  public :: ATMOS_PHY_CH_rn222_tendency

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, private, parameter :: QA_CH = 1

  integer,                public :: ATMOS_PHY_CH_rn222_ntracers = QA_CH

  character(len=H_SHORT), public :: ATMOS_PHY_CH_rn222_NAME(QA_CH) = (/ "RN222" /)
  character(len=H_MID)  , public :: ATMOS_PHY_CH_rn222_DESC(QA_CH) = (/ "Ratio of Rn222 to total mass" /)
  character(len=H_SHORT), public :: ATMOS_PHY_CH_rn222_UNIT(QA_CH) = (/ "Bq/kg" /)

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,  private, parameter :: I_ch_rn222 = 1

  real(RP), private :: ATMOS_PHY_CH_Rn222_decay_ratio ! Decay constant [/s]

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_CH_rn222_setup
    use scale_prc, only: &
       PRC_abort
    implicit none

    real(RP) :: ATMOS_PHY_CH_Rn222_half_life = 3.30048E+5_RP

    namelist / PARAM_ATMOS_PHY_CH_RN222 / &
       ATMOS_PHY_CH_Rn222_half_life

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_CH_rn222_setup",*) 'Setup'
    LOG_INFO("ATMOS_PHY_CH_rn222_setup",*) 'rn222 process'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_CH_RN222,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_PHY_CH_rn222_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_CH_rn222_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_CH_RN222. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_CH_RN222)

    ATMOS_PHY_CH_Rn222_decay_ratio = log(2.0_RP) / ATMOS_PHY_CH_Rn222_half_life

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_CH_rn222_setup",*) 'Characteristics of Rn222'
    LOG_INFO_CONT('(A,E16.6)')             'Half life   [s]   : ', ATMOS_PHY_CH_rn222_half_life
    LOG_INFO_CONT('(A,E16.6)')             'Decay ratio [1/s] : ', ATMOS_PHY_CH_Rn222_decay_ratio

    return
  end subroutine ATMOS_PHY_CH_rn222_setup

  !-----------------------------------------------------------------------------
  !> Chemistry Microphysics
  subroutine ATMOS_PHY_CH_rn222_tendency( &
       KA, KS, KE, &
       IA, IS, IE, &
       JA, JS, JE, &
       QA_CH,      &
       DENS,       &
       QTRC,       &
       RHOQ_t      )
    implicit none

    integer,  intent(in)    :: KA, KS, KE
    integer,  intent(in)    :: IA, IS, IE
    integer,  intent(in)    :: JA, JS, JE
    integer,  intent(in)    :: QA_CH
    real(RP), intent(in)    :: DENS  (KA,IA,JA)
    real(RP), intent(in)    :: QTRC  (KA,IA,JA,QA_CH)
    real(RP), intent(inout) :: RHOQ_t(KA,IA,JA,QA_CH)

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    LOG_PROGRESS(*) 'atmosphere / physics / chemistry / Rn222'

    !--- Decay based on half life
    do j  = JS, JE
    do i  = IS, IE
    do k  = KS, KE
       RHOQ_t(k,i,j,I_ch_rn222) = RHOQ_t(k,i,j,I_ch_rn222) &
                                - DENS(k,i,j) * QTRC(k,i,j,1) * ATMOS_PHY_CH_Rn222_decay_ratio ! [Bq/m3/s]
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_CH_rn222_tendency

end module scale_atmos_phy_ch_rn222
