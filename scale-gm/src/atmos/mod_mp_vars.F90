!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics
!!
!! @par Description
!!         Cloud Microphysics variables
!!
!! @author Team SCALE
!!
!<
! 
!-------------------------------------------------------------------------------
module mod_mp_vars
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_icoA_index

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: mp_vars_setup
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
!  real(RP), public, allocatable :: momz_rmgrid      (:,:,:)
!  real(RP), public, allocatable :: momx_rmgrid      (:,:,:)
!  real(RP), public, allocatable :: momy_rmgrid      (:,:,:)
!  real(RP), public, allocatable :: rhot_rmgrid      (:,:,:)
!  real(RP), public, allocatable :: q_rmgrid         (:,:,:,:)
!  real(RP), public, allocatable :: rho_rmgrid       (:,:,:)
  real(RP), public, allocatable :: CCN_rmgrid       (:,:,:)
  real(RP), public, allocatable :: sflx_rain_rmgrid (:,:)
  real(RP), public, allocatable :: sflx_snow_rmgrid (:,:)
  real(RP), public, allocatable :: Evaporate_rmgrid (:,:,:)

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
  !> Allocate MP variables
  !-----------------------------------------------------------------------------
  subroutine mp_vars_setup
    implicit none

    allocate( CCN_rmgrid       (KA,IA,JA)   )
    allocate( sflx_rain_rmgrid (IA,JA)      )
    allocate( sflx_snow_rmgrid (IA,JA)      )
    allocate( Evaporate_rmgrid (KA,IA,JA)   )

    return
  end subroutine mp_vars_setup

end module mod_mp_vars
