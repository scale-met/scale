!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Radiation
!!
!! @par Description
!!          Common module for Radiation
!!          Heating Rate
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-02-06 (H.Yashiro) [new]
!!
!<
!-------------------------------------------------------------------------------
module scale_atmos_phy_rd_common
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
  public :: ATMOS_PHY_RD_heating

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  ! for direction
  integer,  public, parameter :: I_up = 1
  integer,  public, parameter :: I_dn = 2
  ! for band region
  integer,  public, parameter :: I_LW = 1
  integer,  public, parameter :: I_SW = 2

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
  !> Calc heating rate
  subroutine ATMOS_PHY_RD_heating( &
       flux_rad, &
       DENS,     &
       RHOT,     &
       QTRC,     &
       FZ,       &
       dt,       &
       TEMP_t,   &
       RHOT_t    )
    use scale_atmos_thermodyn, only: &
       THERMODYN_qd   => ATMOS_THERMODYN_qd,   &
       THERMODYN_cv   => ATMOS_THERMODYN_cv,   &
       THERMODYN_rhoe => ATMOS_THERMODYN_rhoe, &
       THERMODYN_rhot => ATMOS_THERMODYN_rhot
    implicit none

    real(RP), intent(in)  :: flux_rad(KA,IA,JA,2,2)
    real(RP), intent(in)  :: DENS    (KA,IA,JA)
    real(RP), intent(in)  :: RHOT    (KA,IA,JA)
    real(RP), intent(in)  :: QTRC    (KA,IA,JA,QA)
    real(RP), intent(in)  :: FZ      (0:KA,IA,JA)
    real(DP), intent(in)  :: dt
    real(RP), intent(out) :: TEMP_t  (KA,IA,JA,3)
    real(RP), intent(out) :: RHOT_t  (KA,IA,JA)

    real(RP) :: RHOE  (KA,IA,JA)
    real(RP) :: RHOE_t(KA,IA,JA,2)
    real(RP) :: RHOT1 (KA,IA,JA)
    real(RP) :: QDRY  (KA,IA,JA)
    real(RP) :: CVtot (KA,IA,JA)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    call THERMODYN_rhoe( RHOE(:,:,:),  & ! [OUT]
                         RHOT(:,:,:),  & ! [IN]
                         QTRC(:,:,:,:) ) ! [IN]

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE

       RHOE_t(k,i,j,I_LW) = ( ( flux_rad(k,i,j,I_LW,I_dn) - flux_rad(k-1,i,j,I_LW,I_dn) ) &
                            - ( flux_rad(k,i,j,I_LW,I_up) - flux_rad(k-1,i,j,I_LW,I_up) ) &
                            ) / ( FZ(k,i,j) - FZ(k-1,i,j) )

       RHOE_t(k,i,j,I_SW) = ( ( flux_rad(k,i,j,I_SW,I_dn) - flux_rad(k-1,i,j,I_SW,I_dn) ) &
                            - ( flux_rad(k,i,j,I_SW,I_up) - flux_rad(k-1,i,j,I_SW,I_up) ) &
                            ) / ( FZ(k,i,j) - FZ(k-1,i,j) )

       RHOE(k,i,j) = RHOE(k,i,j) + dt * ( RHOE_t(k,i,j,I_LW) + RHOE_t(k,i,j,I_SW) )

    enddo
    enddo
    enddo

    call THERMODYN_rhot( RHOT1(:,:,:), & ! [OUT]
                         RHOE(:,:,:),  & ! [IN]
                         QTRC(:,:,:,:) ) ! [IN]

    ! update rhot
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       RHOT_t(k,i,j) = ( RHOT1(k,i,j) - RHOT(k,i,j) ) / dt
    enddo
    enddo
    enddo

    call THERMODYN_qd( QDRY(:,:,:),  & ! [OUT]
                       QTRC(:,:,:,:) ) ! [IN]

    call THERMODYN_cv( CVtot(:,:,:),   & ! [OUT]
                       QTRC (:,:,:,:), & ! [IN]
                       QDRY (:,:,:)    ) ! [IN]

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       TEMP_t(k,i,j,I_LW) = RHOE_t(k,i,j,I_LW) / DENS(k,i,j) / CVtot(k,i,j) * 86400.0_RP ! [K/day]
       TEMP_t(k,i,j,I_SW) = RHOE_t(k,i,j,I_SW) / DENS(k,i,j) / CVtot(k,i,j) * 86400.0_RP ! [K/day]

       TEMP_t(k,i,j,3)    = TEMP_t(k,i,j,I_LW) + TEMP_t(k,i,j,I_SW)
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_RD_heating

end module scale_atmos_phy_rd_common
