!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Radiation
!!
!! @par Description
!!          Atmospheric radiation transfer process wrapper
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-12-06 (S.Nishizawa)   [new]
!!
!<
!-------------------------------------------------------------------------------
module scale_atmos_phy_rd
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
  public :: ATMOS_PHY_RD_setup

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
  abstract interface
     subroutine rd( &
          DENS, RHOT, QTRC,      &
          CZ, FZ,                &
          oceanfrc,              &
          temp_sfc, albedo_land, &
          solins, cosSZA,        &
          flux_rad,              &
          flux_rad_top           )
       use scale_precision
       use scale_grid_index
       use scale_tracer
       implicit none

       real(RP), intent(in)  :: DENS        (KA,IA,JA)
       real(RP), intent(in)  :: RHOT        (KA,IA,JA)
       real(RP), intent(in)  :: QTRC        (KA,IA,JA,QA)
       real(RP), intent(in)  :: CZ          (  KA,IA,JA)    ! UNUSED
       real(RP), intent(in)  :: FZ          (0:KA,IA,JA)
       real(RP), intent(in)  :: oceanfrc    (IA,JA)
       real(RP), intent(in)  :: temp_sfc    (IA,JA)
       real(RP), intent(in)  :: albedo_land (IA,JA,2)
       real(RP), intent(in)  :: solins      (IA,JA)
       real(RP), intent(in)  :: cosSZA      (IA,JA)
       real(RP), intent(out) :: flux_rad    (KA,IA,JA,2,2)
       real(RP), intent(out) :: flux_rad_top(IA,JA,2)
     end subroutine rd
  end interface
  procedure(rd), pointer :: ATMOS_PHY_RD => NULL()
  public :: ATMOS_PHY_RD

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_RD_setup( RD_TYPE )
    use scale_process, only: &
       PRC_MPIstop
#define EXTM(pre, name, post) pre ## name ## post
#define NAME(pre, name, post) EXTM(pre, name, post)
#ifdef RD
    use NAME(scale_atmos_phy_rd_, RD,), only: &
       NAME(ATMOS_PHY_RD_, RD, _setup), &
       NAME(ATMOS_PHY_RD_, RD,)
#else
    use scale_atmos_phy_rd_mstrnx, only: &
       ATMOS_PHY_RD_mstrnx_setup, &
       ATMOS_PHY_RD_mstrnx
#endif
    implicit none

    character(len=H_SHORT), intent(in) :: RD_TYPE
    !---------------------------------------------------------------------------

    select case ( RD_TYPE )
    case ( 'MSTRNX' )
       call ATMOS_PHY_RD_mstrnx_setup( RD_TYPE )
       ATMOS_PHY_RD => ATMOS_PHY_RD_mstrnx
    case default
       write(*,*) 'xxx invalid Radiation type(', trim(RD_TYPE), '). CHECK!'
       call PRC_MPIstop
    end select

    return
  end subroutine ATMOS_PHY_RD_setup

end module scale_atmos_phy_rd
