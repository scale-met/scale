!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Aerosol Microphysics
!!
!! @par Description
!!          Aerosol Microphysics process
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-12-06 (S.Nishizawa)  [new]
!!
!<
!-------------------------------------------------------------------------------
module scale_atmos_phy_ae
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
     subroutine ae( &
          QQA,  &
          DENS, &
          MOMZ, &
          MOMX, &
          MOMY, &
          RHOT, &
          EMIT, &
          NREG, &
          QTRC, &
          CN,   &
          CCN,  &
          RHOQ_t_AE )
       use scale_precision
       use scale_grid_index
       use scale_tracer
       integer,  intent(in)    :: QQA
       real(RP), intent(inout) :: DENS(KA,IA,JA)
       real(RP), intent(inout) :: MOMZ(KA,IA,JA)
       real(RP), intent(inout) :: MOMX(KA,IA,JA)
       real(RP), intent(inout) :: MOMY(KA,IA,JA)
       real(RP), intent(inout) :: RHOT(KA,IA,JA)
       real(RP), intent(inout) :: EMIT(KA,IA,JA,QQA)
       real(RP), intent(in)    :: NREG(KA,IA,JA)
       real(RP), intent(inout) :: QTRC(KA,IA,JA,QA)
       real(RP), intent(out)   :: CN(KA,IA,JA)
       real(RP), intent(out)   :: CCN(KA,IA,JA)
       real(RP), intent(inout) :: RHOQ_t_AE(KA,IA,JA,QA)
     end subroutine ae
     subroutine su
     end subroutine su
     subroutine er( &
          RE, &
          QTRC, &
          RH )
       use scale_precision
       use scale_grid_index
       use scale_tracer
       use scale_atmos_aerosol, only: N_AE
       real(RP), intent(out) :: Re  (KA,IA,JA,N_AE) ! effective radius
       real(RP), intent(in)  :: QTRC(KA,IA,JA,QA)   ! tracer mass concentration [kg/kg]
       real(RP), intent(in)  :: RH  (KA,IA,JA)      ! relative humidity         [0-1]
     end subroutine er
  end interface

  procedure(ae), pointer :: ATMOS_PHY_AE => NULL()
  procedure(su), pointer :: ATMOS_PHY_AE_setup => NULL()
  procedure(er), pointer :: ATMOS_PHY_AE_EffectiveRadius => NULL()
  public :: ATMOS_PHY_AE_config
  public :: ATMOS_PHY_AE
  public :: ATMOS_PHY_AE_setup
  public :: ATMOS_PHY_AE_EffectiveRadius

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public    :: QA_AE  ! number of aerosol microphysical tracers
  integer, public    :: QS_AE  ! start index in QTRC
  integer, public    :: QE_AE  ! end   index in QTRC

  character(len=H_SHORT), pointer, public :: ATMOS_PHY_AE_NAME(:)
  character(len=H_MID),   pointer, public :: ATMOS_PHY_AE_DESC(:)
  character(len=H_SHORT), pointer, public :: ATMOS_PHY_AE_UNIT(:)
  real(RP), pointer, public :: ATMOS_PHY_AE_DENS(:) ! aerosol density [kg/m3]=[g/L]
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
  subroutine ATMOS_PHY_AE_config( AE_TYPE )
    use scale_process, only: &
       PRC_MPIstop
#define EXTM(pre, name, post) pre ## name ## post
#define NAME(pre, name, post) EXTM(pre, name, post)
#ifdef AE
    use NAME(scale_atmos_phy_ae_, AE,), only: &
       NAME(ATMOS_PHY_AE_, AE, _config), &
       NAME(ATMOS_PHY_AE_, AE, _setup), &
       NAME(ATMOS_PHY_AE_, AE,)
       NAME(ATMOS_PHY_AE_, AE, _EffectiveRadius)
#else
    use scale_atmos_phy_ae_dummy, only: &
       ATMOS_PHY_AE_dummy_config, &
       ATMOS_PHY_AE_dummy_setup, &
       ATMOS_PHY_AE_dummy, &
       ATMOS_PHY_AE_dummy_EffectiveRadius, &
       ATMOS_PHY_AE_dummy_NAME, &
       ATMOS_PHY_AE_dummy_UNIT, &
       ATMOS_PHY_AE_dummy_DESC, &
       ATMOS_PHY_AE_dummy_DENS
    use scale_atmos_phy_ae_kajino13, only: &
       ATMOS_PHY_AE_kajino13_config, &
       ATMOS_PHY_AE_kajino13_setup, &
       ATMOS_PHY_AE_kajino13, &
       ATMOS_PHY_AE_kajino13_EffectiveRadius, &
       ATMOS_PHY_AE_kajino13_NAME, &
       ATMOS_PHY_AE_kajino13_UNIT, &
       ATMOS_PHY_AE_kajino13_DESC, &
       ATMOS_PHY_AE_kajino13_DENS
#endif
    implicit none

    character(len=*), intent(in) :: AE_TYPE
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** => ', trim(AE_TYPE), ' is selected.'

    select case( AE_TYPE )
    case ( 'DUMMY', 'NONE' )
       call ATMOS_PHY_AE_dummy_config( &
            AE_TYPE, &
            QA_AE, QS_AE ) ! (out)
       ATMOS_PHY_AE_setup           => ATMOS_PHY_AE_dummy_setup
       ATMOS_PHY_AE                 => ATMOS_PHY_AE_dummy
       ATMOS_PHY_AE_EffectiveRadius => ATMOS_PHY_AE_dummy_EffectiveRadius
       ATMOS_PHY_AE_NAME            => ATMOS_PHY_AE_dummy_NAME
       ATMOS_PHY_AE_DESC            => ATMOS_PHY_AE_dummy_DESC
       ATMOS_PHY_AE_UNIT            => ATMOS_PHY_AE_dummy_UNIT
       ATMOS_PHY_AE_DENS            => ATMOS_PHY_AE_dummy_DENS
    case ( 'KAJINO13' )
       call ATMOS_PHY_AE_kajino13_config( &
            AE_TYPE, &
            QA_AE, QS_AE ) ! (out)
       ATMOS_PHY_AE_setup           => ATMOS_PHY_AE_kajino13_setup
       ATMOS_PHY_AE                 => ATMOS_PHY_AE_kajino13
       ATMOS_PHY_AE_EffectiveRadius => ATMOS_PHY_AE_kajino13_EffectiveRadius
       ATMOS_PHY_AE_NAME            => ATMOS_PHY_AE_kajino13_NAME
       ATMOS_PHY_AE_DESC            => ATMOS_PHY_AE_kajino13_DESC
       ATMOS_PHY_AE_UNIT            => ATMOS_PHY_AE_kajino13_UNIT
       ATMOS_PHY_AE_DENS            => ATMOS_PHY_AE_kajino13_DENS
       write(*,*) '### aerosol type(', AE_TYPE, '). is not recommended in current version!'
    case default
       write(*,*) 'xxx invalid aerosol type(', AE_TYPE, '). CHECK!'
       call PRC_MPIstop
    end select

    QE_AE = QS_AE + QA_AE - 1

    return
  end subroutine ATMOS_PHY_AE_config

end module scale_atmos_phy_ae
