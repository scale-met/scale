!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Surface fluxes
!!
!! @par Description
!!          Flux from/to bottom wall of atmosphere (surface)
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-12-05 (S.Nishizawa)  [new]
!<
!-------------------------------------------------------------------------------
module mod_atmos_phy_sf_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_precision
  use mod_index
  use mod_stdio
  use mod_prof
  use mod_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_SF_driver_setup
  public :: ATMOS_PHY_SF_driver
  public :: ATMOS_PHY_SF_CPL

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

  ! surface flux
  real(RP), allocatable :: SFLX_MOMZ(:,:)
  real(RP), allocatable :: SFLX_MOMX(:,:)
  real(RP), allocatable :: SFLX_MOMY(:,:)
  real(RP), allocatable :: SFLX_POTT(:,:)
  real(RP), allocatable :: SFLX_QV  (:,:)

contains

  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_SF_driver_setup( SF_TYPE )
    use mod_stdio, only: &
       IO_FID_LOG, &
       IO_L, &
       IO_SYSCHR
    use mod_process, only: &
       PRC_MPIstop
    use mod_atmos_phy_sf, only: &
       ATMOS_PHY_SF_setup
    use mod_cpl_vars, only: &
       sw_AtmLnd => CPL_sw_AtmLnd
    implicit none
    character(len=IO_SYSCHR), intent(in) :: SF_TYPE
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[PHY_SURFACEFLUX]/Categ[ATMOS]'

    allocate( SFLX_MOMZ(IA,JA) )
    allocate( SFLX_MOMX(IA,JA) )
    allocate( SFLX_MOMY(IA,JA) )
    allocate( SFLX_POTT(IA,JA) )
    allocate( SFLX_QV  (IA,JA) )


    ! tentative process
    ! finally, surface processes will be located under the coupler.
    if( SF_TYPE /= 'COUPLE' ) then
       call ATMOS_PHY_SF_setup( SF_TYPE )
       call ATMOS_PHY_SF_driver( .true., .false. )
    end if

    return
  end subroutine ATMOS_PHY_SF_driver_setup

  !-----------------------------------------------------------------------------
  ! calculation flux
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_SF_driver( &
       update_flag, &
       history_flag &
       )
    use mod_time, only: &
       dtsf => TIME_DTSEC_ATMOS_PHY_SF, &
       NOWSEC => TIME_NOWDAYSEC
    use mod_const, only: &
       CPdry  => CONST_CPdry,  &
       LH0    => CONST_LH0
    use mod_history, only: &
       HIST_in
    use mod_grid, only: &
       RCDZ => GRID_RCDZ, &
       RFDZ => GRID_RFDZ, &
       CZ   => GRID_CZ
    use mod_atmos_phy_sf, only: &
       ATMOS_PHY_SF
    use mod_atmos_vars, only: &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       RHOT, &
       QTRC, &
       DENS_tp, &
       MOMZ_tp, &
       MOMX_tp, &
       MOMY_tp, &
       RHOT_tp, &
       QTRC_tp
    use mod_ocean_vars, only: &
       SST
    implicit none

    logical, intent(in) :: update_flag
    logical, intent(in), optional :: history_flag

    ! monitor
    real(RP) :: SHFLX(IA,JA) ! sensible heat flux [W/m2]
    real(RP) :: LHFLX(IA,JA) ! latent   heat flux [W/m2]

    integer :: i, j

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Surface flux'

    if ( update_flag ) then
       call ATMOS_PHY_SF( &
            SFLX_MOMZ, SFLX_MOMX, SFLX_MOMY, SFLX_POTT, SFLX_QV, & ! (out)
            DENS, MOMZ, MOMX, MOMY, RHOT, QTRC, SST,             & ! (in)
            CZ, NOWSEC                                           ) ! (in)

       !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
          SHFLX(i,j) = SFLX_POTT(i,j) * CPdry
          LHFLX(i,j) = SFLX_QV  (i,j) * LH0
       end do
       end do

       if ( present(history_flag) ) then
       if ( history_flag ) then
          call HIST_in( SHFLX(:,:), 'SHFLX', 'sensible heat flux', 'W/m2', dtsf )
          call HIST_in( LHFLX(:,:), 'LHFLX', 'latent heat flux',   'W/m2', dtsf )
       end if
       end if

    end if

    !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
       RHOT_tp(KS,i,j) = RHOT_tp(KS,i,j) &
            + ( SFLX_POTT(i,j) &
              + SFLX_QV(i,j) * RHOT(KS,i,j) / DENS(KS,i,j) &
              ) * RCDZ(KS)
       DENS_tp(KS,i,j) = DENS_tp(KS,i,j) &
            + SFLX_QV(i,j) * RCDZ(KS)
       MOMZ_tp(KS,i,j) = MOMZ_tp(KS,i,j) &
            + SFLX_MOMZ(i,j) * RFDZ(KS)
       MOMX_tp(KS,i,j) = MOMX_tp(KS,i,j) &
            + SFLX_MOMX(i,j) * RCDZ(KS)
       MOMY_tp(KS,i,j) = MOMY_tp(KS,i,j) &
            + SFLX_MOMY(i,j) * RCDZ(KS)
       QTRC_tp(KS,i,j,I_QV) = QTRC_tp(KS,i,j,I_QV) &
            + SFLX_QV(i,j) * RCDZ(KS)
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_SF_driver

  subroutine ATMOS_PHY_SF_CPL
    use mod_const, only: &
       CPdry  => CONST_CPdry,  &
       LH0    => CONST_LH0
    use mod_grid, only: &
       RCDZ => GRID_RCDZ, &
       RFDZ => GRID_RFDZ
    use mod_atmos_vars, only: &
       DENS,    &
       RHOT,    &
       DENS_tp, &
       MOMZ_tp, &
       MOMX_tp, &
       MOMY_tp, &
       RHOT_tp, &
       QTRC_tp
    use mod_atmos_vars_sf, only: &
!       SFLX_MOMZ, &
!       SFLX_MOMX, &
!       SFLX_MOMY, &
!       SFLX_SH,   &
!       SFLX_LH,   &
       SFLX_QVAtm
    implicit none

    integer :: i, j

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Surface'

    do j = JS, JE
    do i = IS, IE
!       RHOT_tp(KS,i,j) = RHOT_tp(KS,i,j) &
!            + ( SFLX_SH(i,j)/CPdry &
!              + SFLX_QVAtm(i,j) * RHOT(KS,i,j) / DENS(KS,i,j) &
!              ) * RCDZ(KS)
!       DENS_tp(KS,i,j) = DENS_tp(KS,i,j) &
!            + SFLX_QVAtm(i,j)  * RCDZ(KS)
!       MOMZ_tp(KS,i,j) = MOMZ_tp(KS,i,j) &
!            + SFLX_MOMZ(i,j)   * RFDZ(KS)
!       MOMX_tp(KS,i,j) = MOMX_tp(KS,i,j) &
!            + SFLX_MOMX(i,j)   * RCDZ(KS)
!       MOMY_tp(KS,i,j) = MOMY_tp(KS,i,j) &
!            + SFLX_MOMY(i,j)   * RCDZ(KS)
!       QTRC_tp(KS,i,j,I_QV) = QTRC_tp(KS,i,j,I_QV) &
!            + SFLX_QVAtm(i,j)  * RCDZ(KS)
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_SF_CPL

end module mod_atmos_phy_sf_driver
