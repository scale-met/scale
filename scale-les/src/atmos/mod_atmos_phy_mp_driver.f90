!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics
!!
!! @par Description
!!          Cloud Microphysics wrapper
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-12-06 (S.Nishizawa) [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_phy_mp_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_precision
  use mod_stdio
  use mod_prof
  use mod_grid_index
  use mod_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_MP_driver_setup
  public :: ATMOS_PHY_MP_driver

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
!  real(RP), private, allocatable :: DENS_t(:,:,:)
!  real(RP), private, allocatable :: MOMZ_t(:,:,:)
!  real(RP), private, allocatable :: MOMX_t(:,:,:)
!  real(RP), private, allocatable :: MOMY_t(:,:,:)
!  real(RP), private, allocatable :: RHOT_t(:,:,:)
!  real(RP), private, allocatable :: QTRC_t(:,:,:,:)
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup Cloud Microphysics
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_driver_setup( MP_TYPE )
    use mod_stdio, only: &
       IO_FID_LOG, &
       IO_L, &
       IO_SYSCHR
    use mod_process, only: &
       PRC_MPIstop
    use mod_atmos_phy_mp, only: &
       ATMOS_PHY_MP_setup
    implicit none
    character(len=IO_SYSCHR), intent(in) :: MP_TYPE
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Cloud Microphisics]/Categ[ATMOS]'

!    allocate( DENS_t(KA,IA,JA) )
!    allocate( MOMZ_t(KA,IA,JA) )
!    allocate( MOMX_t(KA,IA,JA) )
!    allocate( MOMY_t(KA,IA,JA) )
!    allocate( RHOT_t(KA,IA,JA) )
!    allocate( QTRC_t(KA,IA,JA,QA) )

    call ATMOS_PHY_MP_setup( MP_TYPE )

!    call ATMOS_PHY_MP_driver( .true., .false. )

    return
  end subroutine ATMOS_PHY_MP_driver_setup

  !-----------------------------------------------------------------------------
  !> Cloud Microphysics
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_driver( update_flag, history_flag )
    use mod_atmos_phy_mp, only: &
       ATMOS_PHY_MP
    use mod_atmos_vars, only: &
       ATMOS_vars_fillhalo, &
       ATMOS_vars_total, &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       RHOT, &
       QTRC
    implicit none
    logical, intent(in) :: update_flag
    logical, intent(in) :: history_flag

    if ( update_flag ) then
       call ATMOS_PHY_MP( &
            DENS, &
            MOMZ, &
            MOMX, &
            MOMY, &
            RHOT, &
            QTRC )

       call ATMOS_vars_fillhalo

       call ATMOS_vars_total
    end if

    return
  end subroutine ATMOS_PHY_MP_driver

end module mod_atmos_phy_mp_driver
