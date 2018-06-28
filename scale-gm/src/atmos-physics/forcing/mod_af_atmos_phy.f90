!-------------------------------------------------------------------------------
!> Module Atmospheric Physics forcing
!!
!! @par Description
!!          This module contains subroutines for physical forcing terms
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_af_atmos_phy
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_atmos_grid_icoA_index

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: af_atmos_phy

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine af_atmos_phy( &
       ijdim,   &
       lat,     &
       lon,     &
       alt,     &
       alth,    &
       rho,     &
       pre,     &
       tem,     &
       vx,      &
       vy,      &
       vz,      &
       w,       &
       q,       &
       ein,     &
       pre_sfc, &
       fvx,     &
       fvy,     &
       fvz,     &
       fw,      &
       fe,      &
       fq,      &
       frho,    &
       precip,  &
       ix,      &
       iy,      &
       iz,      &
       jx,      &
       jy,      &
       jz,      &
       dt       )
    use mod_af_heldsuarez, only: &
       AF_heldsuarez
    use mod_af_dcmip, only: &
       USE_HeldSuarez
    use mod_runconf, only: &
       TRC_VMAX
    use scale_atmos_grid_icoA_index, only: &
       kdim => ADM_kall
    implicit none

    integer,  intent(in)  :: ijdim
    real(RP), intent(in)  :: lat    (ijdim)
    real(RP), intent(in)  :: lon    (ijdim)
    real(RP), intent(in)  :: alt    (ijdim,kdim)
    real(RP), intent(in)  :: alth   (ijdim,kdim)
    real(RP), intent(in)  :: rho    (ijdim,kdim)
    real(RP), intent(in)  :: pre    (ijdim,kdim)
    real(RP), intent(in)  :: tem    (ijdim,kdim)
    real(RP), intent(in)  :: vx     (ijdim,kdim)
    real(RP), intent(in)  :: vy     (ijdim,kdim)
    real(RP), intent(in)  :: vz     (ijdim,kdim)
    real(RP), intent(in)  :: w      (ijdim,kdim)
    real(RP), intent(in)  :: q      (ijdim,kdim,TRC_VMAX)
    real(RP), intent(in)  :: ein    (ijdim,kdim)
    real(RP), intent(in)  :: pre_sfc(ijdim)
    real(RP), intent(out) :: fvx    (ijdim,kdim)
    real(RP), intent(out) :: fvy    (ijdim,kdim)
    real(RP), intent(out) :: fvz    (ijdim,kdim)
    real(RP), intent(out) :: fw     (ijdim,kdim)
    real(RP), intent(out) :: fe     (ijdim,kdim)
    real(RP), intent(out) :: fq     (ijdim,kdim,TRC_VMAX)
    real(RP), intent(out) :: frho   (ijdim,kdim)
    real(RP), intent(out) :: precip (ijdim)
    real(RP), intent(in)  :: ix     (ijdim)
    real(RP), intent(in)  :: iy     (ijdim)
    real(RP), intent(in)  :: iz     (ijdim)
    real(RP), intent(in)  :: jx     (ijdim)
    real(RP), intent(in)  :: jy     (ijdim)
    real(RP), intent(in)  :: jz     (ijdim)
    real(DP), intent(in)  :: dt

    !---------------------------------------------------------------------------

    if ( USE_HeldSuarez ) then
       fw (:,:)   = 0.0_RP
       fq (:,:,:) = 0.0_RP
       frho(:,:)  = 0.0_RP

       call AF_heldsuarez( ijdim,    & ! [IN]
                           lat(:),   & ! [IN]
                           pre(:,:), & ! [IN]
                           tem(:,:), & ! [IN]
                           vx (:,:), & ! [IN]
                           vy (:,:), & ! [IN]
                           vz (:,:), & ! [IN]
                           fvx(:,:), & ! [OUT]
                           fvy(:,:), & ! [OUT]
                           fvz(:,:), & ! [OUT]
                           fe (:,:)  ) ! [OUT]
    endif

    return
  end subroutine af_atmos_phy

end module mod_af_atmos_phy
