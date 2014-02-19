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

  ! surface flux
  real(RP), allocatable :: SFLX_MOMZ(:,:)
  real(RP), allocatable :: SFLX_MOMX(:,:)
  real(RP), allocatable :: SFLX_MOMY(:,:)
  real(RP), allocatable :: SFLX_POTT(:,:)
  real(RP), allocatable :: SFLX_SWU (:,:)
  real(RP), allocatable :: SFLX_LWU (:,:)
  real(RP), allocatable :: SFLX_SH  (:,:)
  real(RP), allocatable :: SFLX_LH  (:,:)
  real(RP), allocatable :: SFLX_QV  (:,:)

contains

  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_SF_driver_setup( SF_TYPE )
    use mod_process, only: &
       PRC_MPIstop
    use mod_atmos_phy_sf, only: &
       ATMOS_PHY_SF_setup
    implicit none

    character(len=H_SHORT), intent(in) :: SF_TYPE
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[PHY_SURFACEFLUX]/Categ[ATMOS]'

    allocate( SFLX_MOMZ(IA,JA) )
    allocate( SFLX_MOMX(IA,JA) )
    allocate( SFLX_MOMY(IA,JA) )
    allocate( SFLX_POTT(IA,JA) )
    allocate( SFLX_SWU (IA,JA) )
    allocate( SFLX_LWU (IA,JA) )
    allocate( SFLX_SH  (IA,JA) )
    allocate( SFLX_LH  (IA,JA) )
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

  subroutine ATMOS_PHY_SF_driver_first
    use mod_const, only: &
       CPdry => CONST_CPdry, &
       RovCP => CONST_RovCP
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
    use mod_cpl_vars, only: &
       CPL_getCPL2Atm
    implicit none

    ! work
    integer :: i, j

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Surface'

    call CPL_getCPL2Atm( &
       SFLX_MOMX, SFLX_MOMY, SFLX_MOMZ, &
       SFLX_SWU, SFLX_LWU,              &
       SFLX_SH, SFLX_LH, SFLX_QV        )

    do j = JS, JE
    do i = IS, IE
       RHOT_tp(KS,i,j) = RHOT_tp(KS,i,j) &
            + ( SFLX_SH(i,j)/CPdry &
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
  end subroutine ATMOS_PHY_SF_driver_first

  subroutine ATMOS_PHY_SF_driver_final
    use mod_const, only: &
       RovCP => CONST_RovCP
    use mod_atmos_vars, only: &
       DENS,    &
       MOMX,    &
       MOMY,    &
       MOMZ,    &
       RHOT,    &
       QTRC
    use mod_atmos_vars_sf, only: &
       PREC, &
       SWD,  &
       LWD
    use mod_atmos_thermodyn, only: &
       temp_pres => ATMOS_THERMODYN_temp_pres
    use mod_cpl_vars, only: &
       CPL_putAtm
    implicit none

    ! work
    integer :: i, j

    real(RP) :: RHOS(IA,JA) ! air density at the sruface [kg/m3]
    real(RP) :: PRES(IA,JA) ! pressure at the surface [Pa]
    real(RP) :: ATMP(IA,JA) ! air temperature at the surface [K]

    real(RP) :: tem(KA,IA,JA) ! temperature [K]
    real(RP) :: pre(KA,IA,JA) ! pressure [Pa]
    real(RP) :: zs (IA,JA)

    call temp_pres( tem (:,:,:),  & ! (out)
                    pre (:,:,:),  & ! (out)
                    DENS(:,:,:),  & ! (in)
                    RHOT(:,:,:),  & ! (in)
                    QTRC(:,:,:,:) ) ! (in)

    zs(:,:) = 0.0_RP
    call sfcval_estimate( RHOS(:,:),   & ! (out)
                          PRES(:,:),   & ! (out)
                          DENS(:,:,:), & ! (in)
                          pre (:,:,:), & ! (in)
                          zs  (:,:)    ) ! (in)

    do j = 1, JA
    do i = 1, IA
      ATMP(i,j) = tem(KS,i,j) * ( PRES(i,j) / pre(KS,i,j) )**RovCP
    end do
    end do

    call CPL_putAtm( &
       DENS(KS,:,:),      &
       MOMX(KS,:,:),      &
       MOMY(KS,:,:),      &
       MOMZ(KS,:,:),      &
       RHOS(:,:),         &
       PRES(:,:),         &
       ATMP(:,:),         &
       QTRC(KS,:,:,I_QV), &
       PREC(:,:),         &
       SWD(:,:),          &
       LWD(:,:)           )

    return
  end subroutine ATMOS_PHY_SF_driver_final

  subroutine sfcval_estimate( &
      sfc_rho, sfc_pre, & ! (out)
      rho, pre, zs      ) ! (in)
    use mod_const, only: &
      GRAV => CONST_GRAV
    use mod_grid, only: &
      CZ => GRID_CZ
    implicit none

    ! argument
    real(RP), intent(out) :: sfc_rho(IA,JA)    ! density at surface [kg/m3]
    real(RP), intent(out) :: sfc_pre(IA,JA)    ! pressure at surface [Pa]
    real(RP), intent(in)  :: rho    (KA,IA,JA) ! density [kg/m3]
    real(RP), intent(in)  :: pre    (KA,IA,JA) ! pressure [Pa]
    real(RP), intent(in)  :: zs     (IA,JA)    ! surface height [m]

    ! work
    integer :: i, j

    real(RP) :: zz
    real(RP) :: z1, z2, z3
    real(RP) :: p1, p2, p3
    real(RP) :: lag_intpl

    lag_intpl( zz, z1, p1, z2, p2, z3, p3 )                &
      = ( (zz-z2) * (zz-z3) ) / ( (z1-z2) * (z1-z3) ) * p1 &
      + ( (zz-z1) * (zz-z3) ) / ( (z2-z1) * (z2-z3) ) * p2 &
      + ( (zz-z1) * (zz-z2) ) / ( (z3-z1) * (z3-z2) ) * p3

    ! estimate surface density (extrapolation)
    do j = 1, JA
    do i = 1, IA
      sfc_rho(i,j) = lag_intpl( zs(i,j),                 &
                                CZ(KS  ), rho(KS  ,i,j), &
                                CZ(KS+1), rho(KS+1,i,j), &
                                CZ(KS+2), rho(KS+2,i,j)  )
    end do
    end do

    ! estimate surface pressure (hydrostatic balance)
    do j = 1, JA
    do i = 1, IA
      sfc_pre(i,j) = pre(KS,i,j)                             &
                   + 0.5_RP * ( sfc_rho(i,j) + rho(KS,i,j) ) &
                   * GRAV * ( CZ(KS) - zs(i,j) )
    end do
    end do

    return
  end subroutine sfcval_estimate

end module mod_atmos_phy_sf_driver
