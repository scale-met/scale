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
    use mod_atmos_dyn, only: &
       ATMOS_DYN_setup
    use mod_atmos_boundary, only: &
       ATMOS_BOUNDARY_setup
    use mod_atmos_phy_mp, only: &
       ATMOS_PHY_MP_setup
    implicit none
    !---------------------------------------------------------------------------

    call ATMOS_vars_setup

    call ATMOS_vars_restart_read

    call ATMOS_REFSTATE_setup

    call ATMOS_DYN_setup

    call ATMOS_BOUNDARY_setup

    call ATMOS_PHY_MP_setup

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
       TIME_NOWSTEP,                     &
       NOWSEC    => TIME_NOWSEC,         &
       TIME_DTSEC_ATMOS_DYN,             &
       TIME_DTSEC_ATMOS_PHY_TB,          &
       TIME_DTSEC_ATMOS_PHY_MP,          &
       do_dyn    => TIME_DOATMOS_DYN,    &
       do_phy_tb => TIME_DOATMOS_PHY_TB, &
       do_phy_mp => TIME_DOATMOS_PHY_MP, &
       do_phy_rd => TIME_DOATMOS_PHY_RD, &
       TIME_rapstart,                    &
       TIME_rapend
    use mod_grid, only: &
       IA   => GRID_IA,   &
       JA   => GRID_JA,   &
       KA   => GRID_KA,   &
       KMAX => GRID_KMAX, &
       IS   => GRID_IS,   &
       IE   => GRID_IE,   &
       JS   => GRID_JS,   &
       JE   => GRID_JE,   &
       KS   => GRID_KS,   &
       KE   => GRID_KE,   &
       WS   => GRID_WS,   &
       WE   => GRID_WE
    use mod_fileio_h, only: &
       FIO_HMID, &
       FIO_REAL8
    use mod_fileio, only: &
       FIO_output
    use mod_atmos_vars, only: &
       QA        => A_QA,            &
       sw_dyn    => ATMOS_sw_dyn,    &
       sw_phy_tb => ATMOS_sw_phy_tb, &
       sw_phy_mp => ATMOS_sw_phy_mp, &
       sw_phy_rd => ATMOS_sw_phy_rd, &
       ATMOS_vars_get,               &
       ATMOS_vars_getall,            &
       ATMOS_vars_put
    use mod_atmos_dyn, only: &
       ATMOS_DYN
    use mod_atmos_boundary, only: &
       ATMOS_BOUNDARY
    use mod_atmos_phy_sf, only: &
       ATMOS_PHY_SF
    use mod_atmos_phy_tb, only: &
       ATMOS_PHY_TB
    use mod_atmos_phy_mp, only: &
       ATMOS_PHY_MP
!    use mod_atmos_phy_rd, only: &
!       ATMOS_PHY_RD
    implicit none

    ! prognostics
    real(8) :: dens(IA,JA,KA)    ! density [kg/m3]
    real(8) :: momx(IA,JA,KA)    ! momentum(x) [kg/m3 * m/s]
    real(8) :: momy(IA,JA,KA)    ! momentum(y) [kg/m3 * m/s]
    real(8) :: momz(IA,JA,KA)    ! momentum(z) [kg/m3 * m/s]
    real(8) :: pott(IA,JA,KA)    ! potential temperature [K]
    real(8) :: qtrc(IA,JA,KA,QA) ! tracer mixing ratio [kg/kg],[1/m3]

    ! diagnostics
    real(8) :: pres(IA,JA,KA)    ! pressure [Pa]
    real(8) :: velx(IA,JA,KA)    ! velocity(x) [m/s]
    real(8) :: vely(IA,JA,KA)    ! velocity(y) [m/s]
    real(8) :: velz(IA,JA,KA)    ! velocity(z) [m/s]
    real(8) :: temp(IA,JA,KA)    ! temperature [K]

    ! surface flux
    real(8) :: FLXij_sfc(IA,JA,3)  ! => FLXij(1:IA,1:JA,WS,1:3,3)
    real(8) :: FLXt_sfc (IA,JA)    ! => FLXt (1:IA,1:JA,WS)
    real(8) :: FLXqv_sfc(IA,JA)    ! => FLXq (1:IA,1:JA,WS,1)

    ! tendency
    real(8) :: dens_t(IA,JA,KA)    ! density [kg/m3]
    real(8) :: momx_t(IA,JA,KA)    ! momentum(x) [kg/m3 * m/s]
    real(8) :: momy_t(IA,JA,KA)    ! momentum(y) [kg/m3 * m/s]
    real(8) :: momz_t(IA,JA,KA)    ! momentum(z) [kg/m3 * m/s]
    real(8) :: pott_t(IA,JA,KA)    ! potential temperature [K]
    real(8) :: qtrc_t(IA,JA,KA,QA) ! tracer mixing ratio   [kg/kg],[1/m3]

    character(len=FIO_HMID) :: desc
    character(len=8)        :: lname
    integer                 :: step = 0
    !---------------------------------------------------------------------------

    if ( mod(TIME_NOWSTEP,10) == 0 ) then
       step = step + 1
       desc = 'temporal history'
       write(lname,'(A,I4.4)') 'ZDEF', KMAX
    endif

    call TIME_rapstart('VARset')
    call ATMOS_vars_get( dens, momx, momy, momz, pott, qtrc )
    call TIME_rapend  ('VARset')


    !########## Dynamics ##########
    call TIME_rapstart('Dynamics')
    if ( sw_dyn .AND. do_dyn ) then 
       call ATMOS_DYN( dens,   momx,   momy,   momz,   pott,   qtrc,  & ! [IN]
                       dens_t, momx_t, momy_t, momz_t, pott_t, qtrc_t ) ! [OUT]

       call ATMOS_BOUNDARY( dens, momx,   momy,   momz,   pott,  & ! [IN]
                                  momx_t, momy_t, momz_t, pott_t ) ! [INOUT]

       dens(:,:,:)   = dens(:,:,:)   + TIME_DTSEC_ATMOS_DYN * dens_t(:,:,:)
       momx(:,:,:)   = momx(:,:,:)   + TIME_DTSEC_ATMOS_DYN * momx_t(:,:,:)
       momy(:,:,:)   = momy(:,:,:)   + TIME_DTSEC_ATMOS_DYN * momy_t(:,:,:)
       momz(:,:,:)   = momz(:,:,:)   + TIME_DTSEC_ATMOS_DYN * momz_t(:,:,:)
       pott(:,:,:)   = pott(:,:,:)   + TIME_DTSEC_ATMOS_DYN * pott_t(:,:,:)
       qtrc(:,:,:,:) = qtrc(:,:,:,:) + TIME_DTSEC_ATMOS_DYN * qtrc_t(:,:,:,:)

       if ( mod(TIME_NOWSTEP,10) == 0 ) then
          call FIO_output( momx_t(IS:IE,JS:JE,KS:KE), 'history', desc, '', 'MOMX_t_dyn', '', '', '', &
                           FIO_REAL8, lname, 1, KMAX, step, NOWSEC, NOWSEC                           )
          call FIO_output( momy_t(IS:IE,JS:JE,KS:KE), 'history', desc, '', 'MOMY_t_dyn', '', '', '', &
                           FIO_REAL8, lname, 1, KMAX, step, NOWSEC, NOWSEC                           )
          call FIO_output( momz_t(IS:IE,JS:JE,KS:KE), 'history', desc, '', 'MOMZ_t_dyn', '', '', '', &
                           FIO_REAL8, lname, 1, KMAX, step, NOWSEC, NOWSEC                           )
          call FIO_output( pott_t(IS:IE,JS:JE,KS:KE), 'history', desc, '', 'POTT_t_dyn', '', '', '', &
                           FIO_REAL8, lname, 1, KMAX, step, NOWSEC, NOWSEC                           )
          call FIO_output( qtrc_t(IS:IE,JS:JE,KS:KE,1), 'history', desc, '', 'QV_t_dyn', '', '', '', &
                           FIO_REAL8, lname, 1, KMAX, step, NOWSEC, NOWSEC                           )
       endif
    endif
    call TIME_rapend  ('Dynamics')

    call TIME_rapstart('VARset')
    call ATMOS_vars_put   ( dens, momx, momy, momz, pott, qtrc  ) ! [IN]
    call ATMOS_vars_getall( dens, momx, momy, momz, pott, qtrc, & ! [OUT]
                            pres, velx, vely, velz, temp        ) ! [OUT]
    call TIME_rapend  ('VARset')


    !########## Turbulence ##########
    call TIME_rapstart('Turbulence')
    if ( sw_phy_tb .AND. do_phy_tb ) then
       call ATMOS_PHY_SF( dens,   pott,   qtrc,          & ! [IN]
                          pres,   velx,   vely,   velz,  & ! [IN]
                          FLXij_sfc, FLXt_sfc, FLXqv_sfc ) ! [OUT]

       call ATMOS_PHY_TB( dens,   pott,   qtrc,                  & ! [IN]
                          velx,   vely,   velz,                  & ! [IN]
                          FLXij_sfc, FLXt_sfc, FLXqv_sfc,        & ! [IN]
                          momx_t, momy_t, momz_t, pott_t, qtrc_t ) ! [OUT]

       momx(:,:,:)   = momx(:,:,:)   + TIME_DTSEC_ATMOS_PHY_TB * momx_t(:,:,:)
       momy(:,:,:)   = momy(:,:,:)   + TIME_DTSEC_ATMOS_PHY_TB * momy_t(:,:,:)
       momz(:,:,:)   = momz(:,:,:)   + TIME_DTSEC_ATMOS_PHY_TB * momz_t(:,:,:)
       pott(:,:,:)   = pott(:,:,:)   + TIME_DTSEC_ATMOS_PHY_TB * pott_t(:,:,:)
       qtrc(:,:,:,:) = qtrc(:,:,:,:) + TIME_DTSEC_ATMOS_PHY_TB * qtrc_t(:,:,:,:)

       if ( mod(TIME_NOWSTEP,10) == 0 ) then
          call FIO_output( momx_t(IS:IE,JS:JE,KS:KE), 'history', desc, '', 'MOMX_t_tb', '', '', '', &
                           FIO_REAL8, lname, 1, KMAX, step, NOWSEC, NOWSEC                          )
          call FIO_output( momy_t(IS:IE,JS:JE,KS:KE), 'history', desc, '', 'MOMY_t_tb', '', '', '', &
                           FIO_REAL8, lname, 1, KMAX, step, NOWSEC, NOWSEC                          )
          call FIO_output( momz_t(IS:IE,JS:JE,KS:KE), 'history', desc, '', 'MOMZ_t_tb', '', '', '', &
                           FIO_REAL8, lname, 1, KMAX, step, NOWSEC, NOWSEC                          )
          call FIO_output( pott_t(IS:IE,JS:JE,KS:KE), 'history', desc, '', 'POTT_t_tb', '', '', '', &
                           FIO_REAL8, lname, 1, KMAX, step, NOWSEC, NOWSEC                          )
          call FIO_output( qtrc_t(IS:IE,JS:JE,KS:KE,1), 'history', desc, '', 'QV_t_tb', '', '', '', &
                           FIO_REAL8, lname, 1, KMAX, step, NOWSEC, NOWSEC                          )
       endif
    endif
    call TIME_rapend  ('Turbulence')


    !########## Microphysics ##########
    call TIME_rapstart('Microphysics')
    if ( sw_phy_mp .AND. do_phy_mp ) then

       call ATMOS_PHY_MP( dens,   momx,   momy,   momz,   pott,   qtrc,  & ! [IN]
                          dens_t, momx_t, momy_t, momz_t, pott_t, qtrc_t ) ! [OUT]

       dens(:,:,:)   = dens(:,:,:)   + TIME_DTSEC_ATMOS_PHY_MP * dens_t(:,:,:)
       momx(:,:,:)   = momx(:,:,:)   + TIME_DTSEC_ATMOS_PHY_MP * momx_t(:,:,:)
       momy(:,:,:)   = momy(:,:,:)   + TIME_DTSEC_ATMOS_PHY_MP * momy_t(:,:,:)
       momz(:,:,:)   = momz(:,:,:)   + TIME_DTSEC_ATMOS_PHY_MP * momz_t(:,:,:)
       pott(:,:,:)   = pott(:,:,:)   + TIME_DTSEC_ATMOS_PHY_MP * pott_t(:,:,:)
       qtrc(:,:,:,:) = qtrc(:,:,:,:) + TIME_DTSEC_ATMOS_PHY_MP * qtrc_t(:,:,:,:)

       if ( mod(TIME_NOWSTEP,10) == 0 ) then
          call FIO_output( dens_t(IS:IE,JS:JE,KS:KE), 'history', desc, '', 'DENS_t_mp', '', '', '', &
                           FIO_REAL8, lname, 1, KMAX, step, NOWSEC, NOWSEC                          )
          call FIO_output( momx_t(IS:IE,JS:JE,KS:KE), 'history', desc, '', 'MOMX_t_mp', '', '', '', &
                           FIO_REAL8, lname, 1, KMAX, step, NOWSEC, NOWSEC                          )
          call FIO_output( momy_t(IS:IE,JS:JE,KS:KE), 'history', desc, '', 'MOMY_t_mp', '', '', '', &
                           FIO_REAL8, lname, 1, KMAX, step, NOWSEC, NOWSEC                          )
          call FIO_output( momz_t(IS:IE,JS:JE,KS:KE), 'history', desc, '', 'MOMZ_t_mp', '', '', '', &
                           FIO_REAL8, lname, 1, KMAX, step, NOWSEC, NOWSEC                          )
          call FIO_output( pott_t(IS:IE,JS:JE,KS:KE), 'history', desc, '', 'POTT_t_mp', '', '', '', &
                           FIO_REAL8, lname, 1, KMAX, step, NOWSEC, NOWSEC                          )
          call FIO_output( qtrc_t(IS:IE,JS:JE,KS:KE,1), 'history', desc, '', 'QV_t_mp', '', '', '', &
                           FIO_REAL8, lname, 1, KMAX, step, NOWSEC, NOWSEC                          )
          call FIO_output( qtrc_t(IS:IE,JS:JE,KS:KE,2), 'history', desc, '', 'QC_t_mp', '', '', '', &
                           FIO_REAL8, lname, 1, KMAX, step, NOWSEC, NOWSEC                          )
          call FIO_output( qtrc_t(IS:IE,JS:JE,KS:KE,3), 'history', desc, '', 'QR_t_mp', '', '', '', &
                           FIO_REAL8, lname, 1, KMAX, step, NOWSEC, NOWSEC                          )
          call FIO_output( qtrc_t(IS:IE,JS:JE,KS:KE,4), 'history', desc, '', 'QI_t_mp', '', '', '', &
                           FIO_REAL8, lname, 1, KMAX, step, NOWSEC, NOWSEC                          )
          call FIO_output( qtrc_t(IS:IE,JS:JE,KS:KE,5), 'history', desc, '', 'QS_t_mp', '', '', '', &
                           FIO_REAL8, lname, 1, KMAX, step, NOWSEC, NOWSEC                          )
          call FIO_output( qtrc_t(IS:IE,JS:JE,KS:KE,6), 'history', desc, '', 'QG_t_mp', '', '', '', &
                           FIO_REAL8, lname, 1, KMAX, step, NOWSEC, NOWSEC                          )
       endif
    endif
    call TIME_rapend  ('Microphysics')


    !########## Radiation ##########
!    call TIME_rapstart('Radiation')
!    if ( sw_phy_rd .AND. do_phy_rd ) then
!       call ATMOS_PHY_RD( dens,   momx,   momy,   momz,   pott,   qtrc,  &
!                          pott_t                                         )
!    endif
!    call TIME_rapend  ('Radiation')


    call TIME_rapstart('VARset')
    call ATMOS_vars_put   ( dens, momx, momy, momz, pott, qtrc  ) ! [IN]
    call ATMOS_vars_getall( dens, momx, momy, momz, pott, qtrc, & ! [OUT]
                            pres, velx, vely, velz, temp        ) ! [OUT]
    call TIME_rapend  ('VARset')

    if ( mod(TIME_NOWSTEP,10) == 0 ) then
       call FIO_output( dens(IS:IE,JS:JE,KS:KE), 'history', desc, '', 'DENS', '', '', '', &
                        FIO_REAL8, lname, 1, KMAX, step, NOWSEC, NOWSEC                   )
       call FIO_output( momx(IS:IE,JS:JE,KS:KE), 'history', desc, '', 'MOMX', '', '', '', &
                        FIO_REAL8, lname, 1, KMAX, step, NOWSEC, NOWSEC                   )
       call FIO_output( momy(IS:IE,JS:JE,KS:KE), 'history', desc, '', 'MOMY', '', '', '', &
                        FIO_REAL8, lname, 1, KMAX, step, NOWSEC, NOWSEC                   )
       call FIO_output( momz(IS:IE,JS:JE,KS:KE), 'history', desc, '', 'MOMZ', '', '', '', &
                        FIO_REAL8, lname, 1, KMAX, step, NOWSEC, NOWSEC                   )
       call FIO_output( pott(IS:IE,JS:JE,KS:KE), 'history', desc, '', 'POTT', '', '', '', &
                        FIO_REAL8, lname, 1, KMAX, step, NOWSEC, NOWSEC                   )
       call FIO_output( pres(IS:IE,JS:JE,KS:KE), 'history', desc, '', 'PRES', '', '', '', &
                        FIO_REAL8, lname, 1, KMAX, step, NOWSEC, NOWSEC                   )
       call FIO_output( velx(IS:IE,JS:JE,KS:KE), 'history', desc, '', 'VELX', '', '', '', &
                        FIO_REAL8, lname, 1, KMAX, step, NOWSEC, NOWSEC                   )
       call FIO_output( vely(IS:IE,JS:JE,KS:KE), 'history', desc, '', 'VELY', '', '', '', &
                        FIO_REAL8, lname, 1, KMAX, step, NOWSEC, NOWSEC                   )
       call FIO_output( velz(IS:IE,JS:JE,KS:KE), 'history', desc, '', 'VELZ', '', '', '', &
                        FIO_REAL8, lname, 1, KMAX, step, NOWSEC, NOWSEC                   )
       call FIO_output( qtrc(IS:IE,JS:JE,KS:KE,1), 'history', desc, '', 'QV', '', '', '', &
                        FIO_REAL8, lname, 1, KMAX, step, NOWSEC, NOWSEC                   )
       call FIO_output( qtrc(IS:IE,JS:JE,KS:KE,2), 'history', desc, '', 'QC', '', '', '', &
                        FIO_REAL8, lname, 1, KMAX, step, NOWSEC, NOWSEC                   )
       call FIO_output( qtrc(IS:IE,JS:JE,KS:KE,3), 'history', desc, '', 'QR', '', '', '', &
                        FIO_REAL8, lname, 1, KMAX, step, NOWSEC, NOWSEC                   )
       call FIO_output( qtrc(IS:IE,JS:JE,KS:KE,4), 'history', desc, '', 'QI', '', '', '', &
                        FIO_REAL8, lname, 1, KMAX, step, NOWSEC, NOWSEC                   )
       call FIO_output( qtrc(IS:IE,JS:JE,KS:KE,5), 'history', desc, '', 'QS', '', '', '', &
                        FIO_REAL8, lname, 1, KMAX, step, NOWSEC, NOWSEC                   )
       call FIO_output( qtrc(IS:IE,JS:JE,KS:KE,6), 'history', desc, '', 'QG', '', '', '', &
                        FIO_REAL8, lname, 1, KMAX, step, NOWSEC, NOWSEC                   )
    endif

    return
  end subroutine ATMOS_step

end module mod_atmos
