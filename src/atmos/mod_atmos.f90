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
!! @li      2011-02-15 (H.Yashiro) [add] Microphysics
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_stdio, only: &
     IO_FID_LOG,  &
     IO_L
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
    use mod_time, only: &
       TIME_DTSEC,                       &
       do_dyn    => TIME_DOATMOS_DYN,    &
       do_phy_tb => TIME_DOATMOS_PHY_TB, &
       do_phy_mp => TIME_DOATMOS_PHY_MP, &
       do_phy_rd => TIME_DOATMOS_PHY_RD, &
       TIME_rapstart,                    &
       TIME_rapend
    use mod_grid, only: &
       KA => GRID_KA, &
       IA => GRID_IA, &
       JA => GRID_JA
    use mod_atmos_vars, only: &
       var => atmos_var,             &
       A_NAME,                       &
       QA  => A_QA,                  &
       I_DENS,                       &
       I_MOMX,                       &
       I_MOMY,                       &
       I_MOMZ,                       &
       I_RHOT,                       &
       I_QC,                         &
       I_QR,                         &
       I_QI,                         &
       I_QS,                         &
       I_QG,                         &
       A_NAME,                       &
       A_DESC,                       &
       A_UNIT,                       &
       sw_dyn    => ATMOS_sw_dyn,    &
       sw_phy_tb => ATMOS_sw_phy_tb, &
       sw_phy_mp => ATMOS_sw_phy_mp, &
       sw_phy_rd => ATMOS_sw_phy_rd, &
       ATMOS_vars_getdiag
    use mod_atmos_dyn, only: &
       ATMOS_DYN
!    use mod_atmos_phy_sf, only: &
!       ATMOS_PHY_SF
!    use mod_atmos_phy_tb, only: &
!       ATMOS_PHY_TB
!    use mod_atmos_phy_mp, only: &
!       ATMOS_PHY_MP
!    use mod_atmos_phy_rd, only: &
!       ATMOS_PHY_RD
    use mod_history, only: &
       HIST_in
    implicit none

    ! diagnostics
    real(8) :: pres(KA,IA,JA)      ! pressure    [Pa]
    real(8) :: velx(KA,IA,JA)      ! velocity(x) [m/s]
    real(8) :: vely(KA,IA,JA)      ! velocity(y) [m/s]
    real(8) :: velz(KA,IA,JA)      ! velocity(z) [m/s]
    real(8) :: temp(KA,IA,JA)      ! temperature [K]
    real(8) :: qtot(KA,IA,JA)      ! Hydrometeor mixing ratio [kg/kg]
    real(8) :: pott(KA,IA,JA)      ! potential temperature [K]

    integer :: iq
    !---------------------------------------------------------------------------

    !########## Dynamics ##########
    call TIME_rapstart('Dynamics')
    if ( sw_dyn .AND. do_dyn ) then
       call ATMOS_DYN
    endif
    call TIME_rapend  ('Dynamics')

    !########## Turbulence ##########

    call TIME_rapstart('Turbulence')
!    if ( sw_phy_tb .AND. do_phy_tb ) then
!       call ATMOS_PHY_SF
!       call ATMOS_PHY_TB
!    endif
    call TIME_rapend  ('Turbulence')

    !########## Microphysics ##########
    call TIME_rapstart('Microphysics')
!    if ( sw_phy_mp .AND. do_phy_mp ) then
!       call ATMOS_PHY_MP
!    endif
    call TIME_rapend  ('Microphysics')

    !########## Radiation ##########
!    call TIME_rapstart('Radiation')
!    if ( sw_phy_rd .AND. do_phy_rd ) then
!       call ATMOS_PHY_RD
!    endif
!    call TIME_rapend  ('Radiation')

    call TIME_rapstart('History')

    call ATMOS_vars_getdiag( pres, velx, vely, velz, temp )

    call HIST_in( var(:,:,:,I_DENS), 'DENS', 'density',     'kg/m3',   '3D', TIME_DTSEC )
    call HIST_in( var(:,:,:,I_MOMZ), 'MOMZ', 'momentum z',  'kg/m2/s', '3D', TIME_DTSEC )
    call HIST_in( var(:,:,:,I_MOMX), 'MOMX', 'momentum x',  'kg/m2/s', '3D', TIME_DTSEC )
    call HIST_in( var(:,:,:,I_MOMY), 'MOMY', 'momentum y',  'kg/m2/s', '3D', TIME_DTSEC )
    call HIST_in( var(:,:,:,I_RHOT), 'RHOT', 'rho * theta', 'kg/m3*K', '3D', TIME_DTSEC )

    if ( QA > 0 ) then
       do iq = 1, QA
          call HIST_in( var(:,:,:,5+iq), A_NAME(5+iq), A_DESC(5+iq), A_UNIT(5+iq), '3D', TIME_DTSEC )
       enddo

       qtot(:,:,:) = var(:,:,:,5+I_QC) &
                   + var(:,:,:,5+I_QR) &
                   + var(:,:,:,5+I_QI) &
                   + var(:,:,:,5+I_QS) &
                   + var(:,:,:,5+I_QG)
       call HIST_in( qtot(:,:,:), 'QTOT', 'Hydrometeor mixing ratio', 'kg/kg', '3D', TIME_DTSEC )
    endif

    call HIST_in( pres(:,:,:), 'PRES', 'pressure',    'Pa',  '3D', TIME_DTSEC )
    call HIST_in( velz(:,:,:), 'W',    'velocity w',  'm/s', '3D', TIME_DTSEC )
    call HIST_in( velx(:,:,:), 'U',    'velocity u',  'm/s', '3D', TIME_DTSEC )
    call HIST_in( vely(:,:,:), 'V',    'velocity v',  'm/s', '3D', TIME_DTSEC )
    call HIST_in( temp(:,:,:), 'T',    'temperature', 'K',   '3D', TIME_DTSEC )

    pott(:,:,:) = var(:,:,:,I_RHOT) / var(:,:,:,I_DENS)
    call HIST_in( pott(:,:,:), 'PT',   'potential temp.', 'K', '3D', TIME_DTSEC )

    call TIME_rapend  ('History')

    return
  end subroutine ATMOS_step

end module mod_atmos
