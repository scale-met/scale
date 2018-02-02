!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Chemistry
!!
!! @par Description
!!          Chemistry process
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2017-02-28 (H.Yashiro)  [new]
!!
!<
!-------------------------------------------------------------------------------
module scale_atmos_phy_ch
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
  abstract interface
     subroutine ch( &
          QQA,   &
          DENS,  &
          QTRC,  &
          RHOQ_t )
       use scale_precision
       use scale_grid_index
       use scale_tracer

       integer,  intent(in)    :: QQA
       real(RP), intent(in)    :: DENS  (KA,IA,JA)
       real(RP), intent(in)    :: QTRC  (KA,IA,JA,QA)
       real(RP), intent(inout) :: RHOQ_t(KA,IA,JA,QQA)
     end subroutine ch

     subroutine su
     end subroutine su
  end interface

  procedure(ch), pointer :: ATMOS_PHY_CH       => NULL()
  procedure(su), pointer :: ATMOS_PHY_CH_setup => NULL()

  public :: ATMOS_PHY_CH_config
  public :: ATMOS_PHY_CH
  public :: ATMOS_PHY_CH_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public :: QA_CH  ! number of chemical tracers
  integer, public :: QS_CH  ! start index in QTRC
  integer, public :: QE_CH  ! end   index in QTRC

  character(len=H_SHORT), pointer, public :: ATMOS_PHY_CH_NAME(:)
  character(len=H_MID),   pointer, public :: ATMOS_PHY_CH_DESC(:)
  character(len=H_SHORT), pointer, public :: ATMOS_PHY_CH_UNIT(:)

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_CH_config( CH_TYPE )
    use scale_process, only: &
       PRC_MPIstop
    use scale_atmos_phy_ch_rn222, only: &
       ATMOS_PHY_CH_rn222_config, &
       ATMOS_PHY_CH_rn222_setup,  &
       ATMOS_PHY_CH_rn222,        &
       ATMOS_PHY_CH_rn222_NAME,   &
       ATMOS_PHY_CH_rn222_UNIT,   &
       ATMOS_PHY_CH_rn222_DESC
    implicit none

    character(len=*), intent(in) :: CH_TYPE
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** => ', trim(CH_TYPE), ' is selected.'

    select case( CH_TYPE )
    case('OFF')
       ! do nothing
    case('RN222')
       call ATMOS_PHY_CH_rn222_config( CH_TYPE, QA_CH, QS_CH )
       ATMOS_PHY_CH_setup => ATMOS_PHY_CH_rn222_setup
       ATMOS_PHY_CH       => ATMOS_PHY_CH_rn222
       ATMOS_PHY_CH_NAME  => ATMOS_PHY_CH_rn222_NAME
       ATMOS_PHY_CH_DESC  => ATMOS_PHY_CH_rn222_DESC
       ATMOS_PHY_CH_UNIT  => ATMOS_PHY_CH_rn222_UNIT
    case default
       write(*,*) 'xxx invalid chemistry type(', CH_TYPE, '). CHECK!'
       call PRC_MPIstop
    end select

    QE_CH = QS_CH + QA_CH - 1

    return
  end subroutine ATMOS_PHY_CH_config

end module scale_atmos_phy_ch
