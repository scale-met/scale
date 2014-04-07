!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Radiation
!!
!! @par Description
!!          dummy code
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-03-21 (H.Yashiro)  [new]
!!
!<
!-------------------------------------------------------------------------------
module scale_atmos_phy_rd_dummy
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
  public :: ATMOS_PHY_RD_dummy_setup
  public :: ATMOS_PHY_RD_dummy

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
  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_RD_dummy_setup( RD_TYPE )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    character(len=H_SHORT), intent(in) :: RD_TYPE
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Physics-RD]/Categ[ATMOS]'
    if( IO_L ) write(IO_FID_LOG,*) '+++ dummy radiation process'

    if ( RD_TYPE /= 'DUMMY' ) then
       if ( IO_L ) write(IO_FID_LOG,*) 'xxx ATMOS_PHY_RD_TYPE is not DUMMY. Check!'
       call PRC_MPIstop
    endif

    return
  end subroutine ATMOS_PHY_RD_dummy_setup

  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_RD_dummy( &
       flux_rad, flux_top, & ! [out]
       solins, cosSZA, & ! [out]
       DENS, RHOT, QTRC, & ! [in]
       CZ, FZ, CDZ, RCDZ, & ! [in]
       REAL_lon, REAL_lat, & ! [in]
       TIME_NOWDATE ) ! [in]
    implicit none
    real(RP), intent(out) :: flux_rad(KA,IA,JA,2,2)
    real(RP), intent(out) :: flux_top(IA,JA,2)
    real(RP), intent(out) :: solins(IA,JA)
    real(RP), intent(out) :: cosSZA(IA,JA)
    real(RP), intent(in)  :: DENS(KA,IA,JA)
    real(RP), intent(in)  :: RHOT(KA,IA,JA)
    real(RP), intent(in)  :: QTRC(KA,IA,JA,QA)
    real(RP), intent(in)  :: CZ(KA)
    real(RP), intent(in)  :: FZ(KA-1)
    real(RP), intent(in)  :: CDZ(KA)
    real(RP), intent(in)  :: RCDZ(KA)
    real(RP), intent(in)  :: REAL_lon(IA,JA)
    real(RP), intent(in)  :: REAL_lat(IA,JA)
    integer , intent(in)  :: TIME_NOWDATE(6)
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Radiation(dummy)'

    return
  end subroutine ATMOS_PHY_RD_dummy

end module scale_atmos_phy_rd_dummy
