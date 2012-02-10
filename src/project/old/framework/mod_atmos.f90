!-------------------------------------------------------------------------------
!> module Atmosphere
!!
!! @par Description
!!          Atmosphere module
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2011-11-11 (H.Yashiro) [new]
!! @li      2011-12-11 (H.Yashiro) [add] Boundary, Surface and Turbulence module
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_setup
  public :: ATMOS_step
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
  !> Setup atmosphere
  !-----------------------------------------------------------------------------
  subroutine ATMOS_setup
    use mod_stdio, only: &
       IO_FID_LOG,  &
       IO_L
    use mod_atmos_vars, only: &
       ATMOS_vars_setup,  &
       ATMOS_vars_restart_read
    use mod_atmos_refstate, only: &
       ATMOS_REFSTATE_setup
    use mod_atmos_boundary, only: &
       ATMOS_BOUNDARY_setup
    use mod_atmos_dyn, only: &
       ATMOS_DYN_setup
!    use mod_atmos_phy_mp, only: &
!       ATMOS_PHY_MP_setup
    implicit none
    !---------------------------------------------------------------------------

    call ATMOS_vars_setup

    call ATMOS_vars_restart_read

    call ATMOS_REFSTATE_setup

    call ATMOS_BOUNDARY_setup

    call ATMOS_DYN_setup

!    call ATMOS_PHY_MP_setup

    return
  end subroutine ATMOS_setup

  !-----------------------------------------------------------------------------
  !> advance atmospheric state
  !-----------------------------------------------------------------------------
  subroutine ATMOS_step
    use mod_stdio, only: &
       IO_FID_LOG,  &
       IO_L
    use mod_time, only: &
       TIME_DTSEC,                       &
       TIME_DTSEC_ATMOS_DYN,             &
       TIME_NSTEP_ATMOS_DYN,             &
       TIME_DTSEC_ATMOS_PHY_TB,          &
       TIME_DTSEC_ATMOS_PHY_MP,          &
       do_dyn    => TIME_DOATMOS_DYN,    &
       do_phy_tb => TIME_DOATMOS_PHY_TB, &
       do_phy_mp => TIME_DOATMOS_PHY_MP, &
       do_phy_rd => TIME_DOATMOS_PHY_RD, &
       TIME_rapstart,                    &
       TIME_rapend
    use mod_grid, only: &
       KA   => GRID_KA, &
       IA   => GRID_IA, &
       JA   => GRID_JA
    use mod_atmos_vars, only: &
       QA        => A_QA,            &
       A_NAME,                       &
       A_DESC,                       &
       A_UNIT,                       &
       sw_dyn    => ATMOS_sw_dyn,    &
       sw_phy_tb => ATMOS_sw_phy_tb, &
       sw_phy_mp => ATMOS_sw_phy_mp, &
       sw_phy_rd => ATMOS_sw_phy_rd, &
       ATMOS_vars_get,               &
       ATMOS_vars_getdiag,           &
       ATMOS_vars_put
    use mod_history, only: &
       HIST_in
    implicit none

    ! prognostics
    real(8) :: dens(KA,IA,JA)      ! density     [kg/m3]
    real(8) :: momx(KA,IA,JA)      ! momentum(x) [kg/m3 * m/s]
    real(8) :: momy(KA,IA,JA)      ! momentum(y) [kg/m3 * m/s]
    real(8) :: momz(KA,IA,JA)      ! momentum(z) [kg/m3 * m/s]
    real(8) :: rhot(KA,IA,JA)      ! rho * theta [kg/m3 * K]
    real(8) :: qtrc(KA,IA,JA,QA)   ! tracer mixing ratio [kg/kg],[1/m3]

    ! diagnostics
    real(8) :: pres(KA,IA,JA)      ! pressure    [Pa]
    real(8) :: velx(KA,IA,JA)      ! velocity(x) [m/s]
    real(8) :: vely(KA,IA,JA)      ! velocity(y) [m/s]
    real(8) :: velz(KA,IA,JA)      ! velocity(z) [m/s]
    real(8) :: temp(KA,IA,JA)      ! temperature [K]
    real(8) :: pott(KA,IA,JA)      ! potential temperature [K]

    ! surface flux
    real(8) :: FLXij_sfc(IA,JA,3)  ! => FLXij(WS,1:IA,1:JA,1:3,3)
    real(8) :: FLXt_sfc (IA,JA)    ! => FLXt (WS,1:IA,1:JA)
    real(8) :: FLXqv_sfc(IA,JA)    ! => FLXq (WS,1:IA,1:JA,I_QV)

    ! tendency
    real(8) :: dens_t(KA,IA,JA)    ! density     [kg/m3]
    real(8) :: momx_t(KA,IA,JA)    ! momentum(x) [kg/m3 * m/s]
    real(8) :: momy_t(KA,IA,JA)    ! momentum(y) [kg/m3 * m/s]
    real(8) :: momz_t(KA,IA,JA)    ! momentum(z) [kg/m3 * m/s]
    real(8) :: rhot_t(KA,IA,JA)    ! rho * theta [kg/m3 * K]
    real(8) :: qtrc_t(KA,IA,JA,QA) ! tracer mixing ratio [kg/kg],[1/m3]

    integer :: iq, step
    !---------------------------------------------------------------------------

    !########## Dynamics ##########
    call TIME_rapstart('Dynamics')
    if ( sw_dyn .AND. do_dyn ) then
       if( IO_L ) write(IO_FID_LOG,*) '+++ Dynamical step'
!       do step = 1, TIME_NSTEP_ATMOS_DYN
!          if( IO_L ) write(IO_FID_LOG,*) '*** Dynamical small step:', step
!       enddo
    endif
    call TIME_rapend  ('Dynamics')

    call TIME_rapstart('Turbulence')
    if ( sw_phy_tb .AND. do_phy_tb ) then
       if( IO_L ) write(IO_FID_LOG,*) '+++ Physical step, TB'
    endif
    call TIME_rapend  ('Turbulence')

    !########## Microphysics ##########
    call TIME_rapstart('Microphysics')
    if ( sw_phy_mp .AND. do_phy_mp ) then
       if( IO_L ) write(IO_FID_LOG,*) '+++ Physical step, MP'
    endif
    call TIME_rapend  ('Microphysics')

    !########## Radiation ##########
    call TIME_rapstart('Radiation')
    if ( sw_phy_rd .AND. do_phy_rd ) then
       if( IO_L ) write(IO_FID_LOG,*) '+++ Physical step, RD'
    endif
    call TIME_rapend  ('Radiation')

    call TIME_rapstart('History')
    call HIST_in( dens(:,:,:), 'DENS', 'density',     'kg/m3',   '3D', TIME_DTSEC )
    call HIST_in( momz(:,:,:), 'MOMZ', 'momentum z',  'kg/m2/s', '3D', TIME_DTSEC )
    call HIST_in( momx(:,:,:), 'MOMX', 'momentum x',  'kg/m2/s', '3D', TIME_DTSEC )
    call HIST_in( momy(:,:,:), 'MOMY', 'momentum y',  'kg/m2/s', '3D', TIME_DTSEC )
    call HIST_in( rhot(:,:,:), 'RHOT', 'rho * theta', 'kg/m3*K', '3D', TIME_DTSEC )
    if ( QA > 0 ) then
       do iq = 1, QA
          call HIST_in( qtrc(:,:,:,iq), A_NAME(5+iq), A_DESC(5+iq), A_UNIT(5+iq), '3D', TIME_DTSEC )
       enddo
    endif

    pott(:,:,:) = rhot(:,:,:) / dens(:,:,:)
    call HIST_in( pott(:,:,:), 'PT',   'potential temp.', 'K', '3D', TIME_DTSEC )

    call HIST_in( pres(:,:,:), 'PRES', 'pressure',    'Pa',  '3D', TIME_DTSEC )
    call HIST_in( velz(:,:,:), 'W',    'velocity w',  'm/s', '3D', TIME_DTSEC )
    call HIST_in( velx(:,:,:), 'U',    'velocity u',  'm/s', '3D', TIME_DTSEC )
    call HIST_in( vely(:,:,:), 'V',    'velocity v',  'm/s', '3D', TIME_DTSEC )
    call HIST_in( temp(:,:,:), 'T',    'temperature', 'K',   '3D', TIME_DTSEC )
    call TIME_rapend  ('History')

    return
  end subroutine ATMOS_step

end module mod_atmos
