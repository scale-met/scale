!-------------------------------------------------------------------------------
!> Program SCALE-LES ver.3
!!
!! @par Description
!!          SCALE: Scalable Computing by Advanced Library and Environment
!!          Numerical model for LES-scale weather
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2011-11-11 (H.Yashiro) [new]
!!
!<
!-------------------------------------------------------------------------------
program scaleles3
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_stdio, only: &
     IO_setup,   &
     IO_FID_LOG, &
     IO_L
  use mod_process, only: &
     PRC_setup,    &
     PRC_MPIstart, &
     PRC_MPIstop
  use mod_const, only: &
     CONST_setup
  use mod_time, only: &
     TIME_setup,           &
     TIME_checkstate,      &
     TIME_advance,         &
     TIME_DOATMOS_step,    &
     TIME_DOOCEAN_step,    &
     TIME_DOATMOS_restart, &
     TIME_DOend,           &
     TIME_rapstart,        &
     TIME_rapend,          &
     TIME_rapreport
  use mod_fileio, only: &
     FIO_setup, &
     FIO_finalize
  use mod_grid, only: &
     GRID_setup
  use mod_comm, only: &
     COMM_setup
  use mod_atmos, only: &
     ATMOS_setup, &
     ATMOS_step
  use mod_atmos_vars, only: &
     ATMOS_vars_restart_write, &
     ATMOS_vars_restart_check, &
     ATMOS_sw_restart,         &
     ATMOS_sw_check
  use mod_ocean, only: &
     OCEAN_setup, &
     OCEAN_step
  use mod_history, only: &
     HIST_setup, &
     HIST_write
  !-----------------------------------------------------------------------------
  implicit none
  !-----------------------------------------------------------------------------
  !
  !++ parameters & variables
  !
  !=============================================================================

  !########## Initial setup ##########

  ! setup standard I/O
  call IO_setup

  ! start MPI
  call PRC_MPIstart

  ! setup process
  call PRC_setup

  ! setup constants
  call CONST_setup

  ! setup time
  call TIME_setup
  call TIME_rapstart('Initialize')

  ! setup file I/O
  call FIO_setup

  ! setup horisontal/veritical grid system
  call GRID_setup

  ! setup mpi communication
  call COMM_setup

  ! setup history
  call HIST_setup

  ! setup atmosphere
  call ATMOS_vars_setup
  call ATMOS_vars_restart_read
  call ATMOS_REFSTATE_setup
  call ATMOS_BOUNDARY_setup
  call ATMOS_DYN_setup

  ! setup ocean
  call OCEAN_vars_setup
  allocate( ocean_var(IA,JA,1,O_VA) )
  ocean_var(:,:,:,:) = OCEAN_SST ! tentative: put contstant value

  call TIME_rapend('Initialize')


  !########## main ##########

  if( IO_L ) write(IO_FID_LOG,*)
  if( IO_L ) write(IO_FID_LOG,*) '++++++ START TIMESTEP ++++++'
  call TIME_rapstart('Main Loop(Total)')

  do

    call TIME_checkstate

    !########## Dynamics ##########
    call TIME_rapstart('Dynamics')
    if ( sw_dyn .AND. do_dyn ) then
       call ATMOS_DYN( &
               var, var_s, diagvar, &
               CDZ, CDX, CDY, RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY, &
               ray_damp, DAMP_var, DAMP_alpha, &
               num_diff, dens_diff, pott_diff, REF_dens, REF_pott, &
               mflx_hi, mflx_lo, qflx_hi, qflx_lo, qflx_anti, ddiv, rjpls, rjmns, &
               KA, IA, JA, VA, QA, &
               KS, KE, WS, WE, IS, IE, JS, JE )
    endif
    call TIME_rapend  ('Dynamics')

    call TIME_rapstart('VARset')
    call ATMOS_vars_get( dens, momx, momy, momz, rhot, qtrc )
    call ATMOS_vars_getdiag( pres, velx, vely, velz, temp )
    call TIME_rapend  ('VARset')

    !########## Turbulence ##########

!    call TIME_rapstart('Turbulence')
!    if ( sw_phy_tb .AND. do_phy_tb ) then
!       momx_t(:,:,:)   = 0.D0
!       momy_t(:,:,:)   = 0.D0
!       momz_t(:,:,:)   = 0.D0
!       rhot_t(:,:,:)   = 0.D0
!       qtrc_t(:,:,:,:) = 0.D0

!       call ATMOS_PHY_SF( dens,   rhot,   qtrc,          & ! [IN]
!                          pres,   velx,   vely,   velz,  & ! [IN]
!                          FLXij_sfc, FLXt_sfc, FLXqv_sfc ) ! [OUT]

!       call ATMOS_PHY_TB( dens,   rhot,   qtrc,                  & ! [IN]
!                          velx,   vely,   velz,                  & ! [IN]
!                          FLXij_sfc, FLXt_sfc, FLXqv_sfc,        & ! [IN]
!                          momx_t, momy_t, momz_t, rhot_t, qtrc_t ) ! [OUT]

!       momx(:,:,:)   = momx(:,:,:)   + TIME_DTSEC_ATMOS_PHY_TB * momx_t(:,:,:)
!       momy(:,:,:)   = momy(:,:,:)   + TIME_DTSEC_ATMOS_PHY_TB * momy_t(:,:,:)
!       momz(:,:,:)   = momz(:,:,:)   + TIME_DTSEC_ATMOS_PHY_TB * momz_t(:,:,:)
!       rhot(:,:,:)   = rhot(:,:,:)   + TIME_DTSEC_ATMOS_PHY_TB * rhot_t(:,:,:)
!       qtrc(:,:,:,:) = qtrc(:,:,:,:) + TIME_DTSEC_ATMOS_PHY_TB * qtrc_t(:,:,:,:)

!    endif
!    call TIME_rapend  ('Turbulence')

    !########## Microphysics ##########
!    call TIME_rapstart('Microphysics')
!    if ( sw_phy_mp .AND. do_phy_mp ) then

!       call ATMOS_PHY_MP( dens,   momx,   momy,   momz,   rhot,   qtrc,  & ! [IN]
!                          dens_t, momx_t, momy_t, momz_t, rhot_t, qtrc_t ) ! [OUT]

!       dens(:,:,:)   = dens(:,:,:)   + TIME_DTSEC_ATMOS_PHY_MP * dens_t(:,:,:)
!       momx(:,:,:)   = momx(:,:,:)   + TIME_DTSEC_ATMOS_PHY_MP * momx_t(:,:,:)
!       momy(:,:,:)   = momy(:,:,:)   + TIME_DTSEC_ATMOS_PHY_MP * momy_t(:,:,:)
!       momz(:,:,:)   = momz(:,:,:)   + TIME_DTSEC_ATMOS_PHY_MP * momz_t(:,:,:)
!       rhot(:,:,:)   = rhot(:,:,:)   + TIME_DTSEC_ATMOS_PHY_MP * rhot_t(:,:,:)
!       qtrc(:,:,:,:) = qtrc(:,:,:,:) + TIME_DTSEC_ATMOS_PHY_MP * qtrc_t(:,:,:,:)

!    endif
!    call TIME_rapend  ('Microphysics')


    !########## Radiation ##########
!    call TIME_rapstart('Radiation')
!    if ( sw_phy_rd .AND. do_phy_rd ) then
!       call ATMOS_PHY_RD( dens,   momx,   momy,   momz,   rhot,   qtrc,  &
!                          rhot_t                                         )
!    endif
!    call TIME_rapend  ('Radiation')


!    call TIME_rapstart('VARset')
!    call ATMOS_vars_put   ( dens, momx, momy, momz, rhot, qtrc  ) ! [IN]
!    call TIME_rapend  ('VARset')

    call TIME_rapstart('History')
    call HIST_in( dens(:,:,:), 'DENS', 'density',     'kg/m3',   '3D', TIME_DTSEC )
    call HIST_in( momz(:,:,:), 'MOMZ', 'momentum z',  'kg/m2/s', '3D', TIME_DTSEC )
    call HIST_in( momx(:,:,:), 'MOMX', 'momentum x',  'kg/m2/s', '3D', TIME_DTSEC )
    call HIST_in( momy(:,:,:), 'MOMY', 'momentum y',  'kg/m2/s', '3D', TIME_DTSEC )
    call HIST_in( rhot(:,:,:), 'RHOT', 'rho * theta', 'kg/m3*K', '3D', TIME_DTSEC )

    pott(:,:,:) = rhot(:,:,:) / dens(:,:,:)
    call HIST_in( pott(:,:,:), 'PT',   'potential temp.', 'K', '3D', TIME_DTSEC )

    call HIST_in( pres(:,:,:), 'PRES', 'pressure',    'Pa',  '3D', TIME_DTSEC )
    call HIST_in( velz(:,:,:), 'W',    'velocity w',  'm/s', '3D', TIME_DTSEC )
    call HIST_in( velx(:,:,:), 'U',    'velocity u',  'm/s', '3D', TIME_DTSEC )
    call HIST_in( vely(:,:,:), 'V',    'velocity v',  'm/s', '3D', TIME_DTSEC )
    call HIST_in( temp(:,:,:), 'T',    'temperature', 'K',   '3D', TIME_DTSEC )
    call TIME_rapend  ('History')

    ! time advance
    call TIME_advance

    ! history file output
    call HIST_write

    ! restart output
    if ( ATMOS_sw_restart .AND. TIME_DOATMOS_restart ) call ATMOS_vars_restart_write

    if ( TIME_DOend ) exit

  enddo

  call TIME_rapend('Main Loop(Total)')
  if( IO_L ) write(IO_FID_LOG,*) '++++++ END TIMESTEP ++++++'
  if( IO_L ) write(IO_FID_LOG,*)


  !########## Finalize ##########

  ! check data
  if ( ATMOS_sw_check ) call ATMOS_vars_restart_check

  call TIME_rapreport

  call FIO_finalize
  ! stop MPI
  call PRC_MPIstop

  stop
  !=============================================================================
end program scaleles3
