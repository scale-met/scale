!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Surface fluxes
!!
!! @par Description
!!          Flux from/to bottom boundary of atmosphere (surface)
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-12-05 (S.Nishizawa)  [new]
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module mod_atmos_phy_sf_driver
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
  public :: ATMOS_PHY_SF_driver_setup
  public :: ATMOS_PHY_SF_driver
  public :: ATMOS_PHY_SF_driver_first
  public :: ATMOS_PHY_SF_driver_final

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
  !> Setup
  subroutine ATMOS_PHY_SF_driver_setup
    use scale_atmos_phy_sf, only: &
       ATMOS_PHY_SF_setup
    use mod_atmos_admin, only: &
       ATMOS_PHY_SF_TYPE, &
       ATMOS_sw_phy_sf
    use mod_atmos_phy_sf_vars, only: &
       SFC_beta => ATMOS_PHY_SF_SFC_beta, &
       SFC_Z0   => ATMOS_PHY_SF_SFC_Z0
    use mod_cpl_vars, only: &
       CPL_sw => CPL_sw_ALL
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[ATMOS PHY_SF] / Origin[SCALE-LES]'

    if ( ATMOS_sw_phy_sf ) then

       ! setup library component
       call ATMOS_PHY_SF_setup( ATMOS_PHY_SF_TYPE )

       if ( .NOT. CPL_sw ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Coupler is disabled.'
          if( IO_L ) write(IO_FID_LOG,*) '*** SFC_beta is assumed to be 1.'
          if( IO_L ) write(IO_FID_LOG,*) '*** SFC_Z0   is assumed to be 0.'
          SFC_beta(:,:) = 1.0_RP
          SFC_Z0  (:,:) = 0.0_RP
       endif

       ! run once (only for the diagnostic value)
       call ATMOS_PHY_SF_driver( .true., .false. )

    else

       if( IO_L ) write(IO_FID_LOG,*) '*** ATMOS_PHY_SF is disabled.'
       if( IO_L ) write(IO_FID_LOG,*) '*** SFC_TEMP, SFC_albedo, SFC_albedo_land is set in ATMOS_PHY_SF_vars.'

    endif

    return
  end subroutine ATMOS_PHY_SF_driver_setup

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine ATMOS_PHY_SF_driver( update_flag, history_flag )
    use scale_const, only: &
       CPdry => CONST_CPdry,  &
       LH0   => CONST_LH0
    use scale_grid, only: &
       RCDZ => GRID_RCDZ, &
       RFDZ => GRID_RFDZ, &
       CZ   => GRID_CZ
    use scale_time, only: &
       dt_SF  => TIME_DTSEC_ATMOS_PHY_SF, &
       NOWSEC => TIME_NOWDAYSEC
    use scale_history, only: &
       HIST_in
    use scale_atmos_phy_sf, only: &
       ATMOS_PHY_SF
    use mod_atmos_vars, only: &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       RHOT, &
       QTRC, &
       DENS_t => DENS_tp, &
       MOMZ_t => MOMZ_tp, &
       MOMX_t => MOMX_tp, &
       MOMY_t => MOMY_tp, &
       RHOT_t => RHOT_tp, &
       QTRC_t => QTRC_tp
    use mod_cpl_vars, only: &
       SST
    use mod_atmos_phy_sf_vars, only: &
       DENS_t_SF => ATMOS_PHY_SF_DENS_t,  &
       MOMZ_t_SF => ATMOS_PHY_SF_MOMZ_t,  &
       MOMX_t_SF => ATMOS_PHY_SF_MOMX_t,  &
       MOMY_t_SF => ATMOS_PHY_SF_MOMY_t,  &
       RHOT_t_SF => ATMOS_PHY_SF_RHOT_t,  &
       QTRC_t_SF => ATMOS_PHY_SF_QTRC_t,  &
       ZMFLX     => ATMOS_PHY_SF_ZMFLX,   &
       XMFLX     => ATMOS_PHY_SF_XMFLX,   &
       YMFLX     => ATMOS_PHY_SF_YMFLX,   &
       POTTFLX   => ATMOS_PHY_SF_POTTFLX, &
       SHFLX     => ATMOS_PHY_SF_SHFLX,   &
       LHFLX     => ATMOS_PHY_SF_LHFLX,   &
       QVFLX     => ATMOS_PHY_SF_QVFLX
    implicit none

    logical, intent(in) :: update_flag
    logical, intent(in) :: history_flag

    integer :: i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Surface flux'

    if ( update_flag ) then
       call ATMOS_PHY_SF( &
            ZMFLX, XMFLX, YMFLX, POTTFLX, QVFLX,     & ! (out)
            DENS, MOMZ, MOMX, MOMY, RHOT, QTRC, SST, & ! (in)
            CZ, NOWSEC                               ) ! (in)

       do j = JS, JE
       do i = IS, IE
          SHFLX(i,j) = POTTFLX(i,j) * CPdry
          LHFLX(i,j) = QVFLX  (i,j) * LH0
       enddo
       enddo

       if ( history_flag ) then
          call HIST_in( SHFLX(:,:), 'SHFLX', 'sensible heat flux', 'W/m2', dt_SF )
          call HIST_in( LHFLX(:,:), 'LHFLX', 'latent heat flux',   'W/m2', dt_SF )
       endif

    endif

    !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
       RHOT_t(KS,i,j) = RHOT_t(KS,i,j) &
            + ( POTTFLX(i,j) &
            + QVFLX(i,j) * RHOT(KS,i,j) / DENS(KS,i,j) &
              ) * RCDZ(KS)
       DENS_t(KS,i,j) = DENS_t(KS,i,j) &
            + QVFLX(i,j) * RCDZ(KS)
       MOMZ_t(KS,i,j) = MOMZ_t(KS,i,j) &
            + ZMFLX(i,j) * RFDZ(KS)
       MOMX_t(KS,i,j) = MOMX_t(KS,i,j) &
            + XMFLX(i,j) * RCDZ(KS)
       MOMY_t(KS,i,j) = MOMY_t(KS,i,j) &
            + YMFLX(i,j) * RCDZ(KS)
       QTRC_t(KS,i,j,I_QV) = QTRC_t(KS,i,j,I_QV) &
            + QVFLX(i,j) * RCDZ(KS)
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_SF_driver

  subroutine ATMOS_PHY_SF_driver_first
    use scale_const, only: &
       CPdry => CONST_CPdry, &
       RovCP => CONST_RovCP
    use scale_grid, only: &
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
    use mod_cpl_vars, only: &
       CPL_getCPL2Atm
    use mod_atmos_phy_sf_vars, only: &
       DENS_t_SF => ATMOS_PHY_SF_DENS_t,  &
       MOMZ_t_SF => ATMOS_PHY_SF_MOMZ_t,  &
       MOMX_t_SF => ATMOS_PHY_SF_MOMX_t,  &
       MOMY_t_SF => ATMOS_PHY_SF_MOMY_t,  &
       RHOT_t_SF => ATMOS_PHY_SF_RHOT_t,  &
       QTRC_t_SF => ATMOS_PHY_SF_QTRC_t,  &
       ZMFLX     => ATMOS_PHY_SF_ZMFLX,   &
       XMFLX     => ATMOS_PHY_SF_XMFLX,   &
       YMFLX     => ATMOS_PHY_SF_YMFLX,   &
       SWUFLX    => ATMOS_PHY_SF_SWUFLX,  &
       LWUFLX    => ATMOS_PHY_SF_LWUFLX,  &
       SHFLX     => ATMOS_PHY_SF_SHFLX,   &
       LHFLX     => ATMOS_PHY_SF_LHFLX,   &
       QVFLX     => ATMOS_PHY_SF_QVFLX
    implicit none

    ! work
    integer :: i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Surface'

    call CPL_getCPL2Atm( &
       XMFLX, YMFLX, ZMFLX, &
       SWUFLX, LWUFLX,      &
       SHFLX, LHFLX, QVFLX  )

    do j = JS, JE
    do i = IS, IE
       RHOT_tp(KS,i,j) = RHOT_tp(KS,i,j) &
            + ( SHFLX(i,j) / CPdry &
            + QVFLX(i,j) * RHOT(KS,i,j) / DENS(KS,i,j) &
              ) * RCDZ(KS)
       DENS_tp(KS,i,j) = DENS_tp(KS,i,j) &
            + QVFLX(i,j) * RCDZ(KS)
       MOMZ_tp(KS,i,j) = MOMZ_tp(KS,i,j) &
            + ZMFLX(i,j) * RFDZ(KS)
       MOMX_tp(KS,i,j) = MOMX_tp(KS,i,j) &
            + XMFLX(i,j) * RCDZ(KS)
       MOMY_tp(KS,i,j) = MOMY_tp(KS,i,j) &
            + YMFLX(i,j) * RCDZ(KS)
       QTRC_tp(KS,i,j,I_QV) = QTRC_tp(KS,i,j,I_QV) &
            + QVFLX(i,j) * RCDZ(KS)
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_SF_driver_first

  subroutine ATMOS_PHY_SF_driver_final
    use scale_const, only: &
       RovCP => CONST_RovCP
    use scale_atmos_thermodyn, only: &
       temp_pres => ATMOS_THERMODYN_temp_pres
    use mod_atmos_vars, only: &
       DENS,    &
       MOMX,    &
       MOMY,    &
       MOMZ,    &
       RHOT,    &
       QTRC
    use mod_cpl_vars, only: &
       CPL_putAtm
    use mod_atmos_phy_sf_vars, only: &
       PREC => ATMOS_PHY_SF_PREC, &
       SWD  => ATMOS_PHY_SF_SWD,  &
       LWD  => ATMOS_PHY_SF_LWD,  &
       ATMOS_PHY_SF_vars_fillhalo
    implicit none

    ! work
    integer :: i, j

    real(RP) :: RHOS(IA,JA) ! air density at the sruface [kg/m3]
    real(RP) :: PRES(IA,JA) ! pressure at the surface [Pa]
    real(RP) :: TMPS(IA,JA) ! air temperature at the surface [K]

    real(RP) :: tem(KA,IA,JA) ! temperature [K]
    real(RP) :: pre(KA,IA,JA) ! pressure [Pa]
    !---------------------------------------------------------------------------

    call temp_pres( tem (:,:,:),  & ! (out)
                    pre (:,:,:),  & ! (out)
                    DENS(:,:,:),  & ! (in)
                    RHOT(:,:,:),  & ! (in)
                    QTRC(:,:,:,:) ) ! (in)

    call sfcval_estimate( RHOS(:,:),   & ! (out)
                          PRES(:,:),   & ! (out)
                          DENS(:,:,:), & ! (in)
                          pre (:,:,:)  ) ! (in)

    do j = 1, JA
    do i = 1, IA
      TMPS(i,j) = tem(KS,i,j) * ( PRES(i,j) / pre(KS,i,j) )**RovCP
    end do
    end do

    call ATMOS_PHY_SF_vars_fillhalo

    call CPL_putAtm( &
       DENS(KS,:,:),      &
       MOMX(KS,:,:),      &
       MOMY(KS,:,:),      &
       MOMZ(KS,:,:),      &
       RHOS(:,:),         &
       PRES(:,:),         &
       TMPS(:,:),         &
       tem(KS,:,:),       &
       QTRC(KS,:,:,I_QV), &
       PREC(:,:),         &
       SWD(:,:),          &
       LWD(:,:)           )

    return
  end subroutine ATMOS_PHY_SF_driver_final

  subroutine sfcval_estimate( &
      sfc_rho, sfc_pre, & ! (out)
      rho, pre          ) ! (in)
    use scale_const, only: &
      GRAV => CONST_GRAV
    use scale_grid_real, only: &
      CZ => REAL_CZ, &
      FZ => REAL_FZ
    implicit none

    ! argument
    real(RP), intent(out) :: sfc_rho(IA,JA)    ! density at surface [kg/m3]
    real(RP), intent(out) :: sfc_pre(IA,JA)    ! pressure at surface [Pa]
    real(RP), intent(in)  :: rho    (KA,IA,JA) ! density [kg/m3]
    real(RP), intent(in)  :: pre    (KA,IA,JA) ! pressure [Pa]

    ! work
    integer :: i, j
    !---------------------------------------------------------------------------

    ! estimate surface density (extrapolation)
    do j = 1, JA
    do i = 1, IA
      sfc_rho(i,j) = lag_intpl( FZ(KS-1,i,j),                &
                                CZ(KS  ,i,j), rho(KS  ,i,j), &
                                CZ(KS+1,i,j), rho(KS+1,i,j), &
                                CZ(KS+2,i,j), rho(KS+2,i,j)  )
    end do
    end do

    ! estimate surface pressure (hydrostatic balance)
    do j = 1, JA
    do i = 1, IA
      sfc_pre(i,j) = pre(KS,i,j)                             &
                   + 0.5_RP * ( sfc_rho(i,j) + rho(KS,i,j) ) &
                   * GRAV * ( CZ(KS,i,j) - FZ(KS-1,i,j) )
    end do
    end do

    return
  end subroutine sfcval_estimate

  function lag_intpl( zz, z1, p1, z2, p2, z3, p3 )
    implicit none

    ! argument
    real(RP), intent(in) :: zz
    real(RP), intent(in) :: z1, z2, z3
    real(RP), intent(in) :: p1, p2, p3
    ! function
    real(RP) :: lag_intpl

    lag_intpl &
      = ( (zz-z2) * (zz-z3) ) / ( (z1-z2) * (z1-z3) ) * p1 &
      + ( (zz-z1) * (zz-z3) ) / ( (z2-z1) * (z2-z3) ) * p2 &
      + ( (zz-z1) * (zz-z2) ) / ( (z3-z1) * (z3-z2) ) * p3

    return
  end function lag_intpl

end module mod_atmos_phy_sf_driver
