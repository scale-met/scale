!-------------------------------------------------------------------------------
!> Program make tool for initial states for SCALE-LES ver.3
!!
!! @par Description
!!          SCALE: Scalable Computing by Advanced Library and Environment
!!          Numerical model for LES-scale weather
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2012-04-08 (H.Yashiro)  [mod] merge all init programs
!!
!<
!-------------------------------------------------------------------------------
program scaleinit
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
  use mod_random, only: &
     RANDOM_setup
  use mod_time, only: &
     TIME_setup,    &
     TIME_rapstart, &
     TIME_rapend,   &
     TIME_rapreport
  use mod_fileio, only: &
     FIO_setup, &
     FIO_finalize
  use mod_grid, only: &
     GRID_setup
  use mod_geometrics, only: &
     GEOMETRICS_setup
  use mod_comm, only: &
     COMM_setup
  use mod_history, only: &
     HIST_setup, &
     HIST_write
  use mod_monitor, only: &
     MONIT_setup, &
     MONIT_write, &
     MONIT_finalize
  use mod_atmos_vars, only: &
     ATMOS_vars_setup, &
     ATMOS_vars_fillhalo, &
     ATMOS_vars_restart_write
  use mod_atmos_thermodyn, only: &
     ATMOS_THERMODYN_setup
  use mod_mkinit, only: &
     MKINIT_TYPE,     &
     I_PLANESTATE,    &
     I_TRACERBUBBLE,  &
     I_COLDBUBBLE,    &
     I_WARMBUBBLE,    &
     I_KHWAVE,        &
     I_TURBULENCE,    &
     I_SUPERCELL,     &
     I_SQUALLINE,     &
     I_DYCOMS2_RF01,  &
     I_DYCOMS2_RF02,  &
     I_DYCOMS2_RF01_hbinw,&
     I_WARMBUBBLE_hbinw,  &
     MKINIT_setup,        &
     MKINIT_planestate,   &
     MKINIT_tracerbubble, &
     MKINIT_coldbubble,   &
     MKINIT_warmbubble,   &
     MKINIT_khwave,       &
     MKINIT_turbulence,   &
     MKINIT_supercell,    &
     MKINIT_squalline,    &
     MKINIT_DYCOMS2_RF01, &
     MKINIT_DYCOMS2_RF02, &
     MKINIT_DYCOMS2_RF01_hbinw, &
     MKINIT_WARMBUBBLE_hbinw
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

  ! setup random number
  call RANDOM_setup

  ! setup time
  call TIME_setup
  call TIME_rapstart('Initialize')

  ! setup file I/O
  call FIO_setup

  ! setup horisontal/veritical grid system
  call GRID_setup

  ! setup geometrics
  call GEOMETRICS_setup

  ! setup mpi communication
  call COMM_setup

  ! setup history
  call HIST_setup

  ! setup monitor
  call MONIT_setup

  call TIME_rapend('Initialize')


  !########## main ##########

  if( IO_L ) write(IO_FID_LOG,*)
  if( IO_L ) write(IO_FID_LOG,*) '++++++ START MAKING INITIAL DATA ++++++'
  call TIME_rapstart('Main')

  ! setup restart
  call ATMOS_THERMODYN_setup
  call ATMOS_vars_setup

  ! setup mkinit
  call MKINIT_setup

  select case(MKINIT_TYPE)
  case(I_PLANESTATE)
     call MKINIT_planestate
  case(I_TRACERBUBBLE)
     call MKINIT_tracerbubble
  case(I_COLDBUBBLE)
     call MKINIT_coldbubble
  case(I_WARMBUBBLE)
     call MKINIT_warmbubble
  case(I_KHWAVE)
     call MKINIT_khwave
  case(I_TURBULENCE)
     call MKINIT_turbulence
  case(I_SUPERCELL)
     call MKINIT_supercell
  case(I_SQUALLINE)
     call MKINIT_squalline
  case(I_DYCOMS2_RF01)
     call MKINIT_DYCOMS2_RF01
  case(I_DYCOMS2_RF02)
     call MKINIT_DYCOMS2_RF02
  case(I_DYCOMS2_RF01_hbinw)
     call MKINIT_DYCOMS2_RF01_hbinw
  case(I_WARMBUBBLE_hbinw)
     call MKINIT_WARMBUBBLE_hbinw
  case default
     write(*,*) ' xxx Unsupported TYPE:', MKINIT_TYPE
     call PRC_MPIstop
  endselect

  call ATMOS_vars_fillhalo
  ! output restart
  call ATMOS_vars_restart_write

  call TIME_rapend('Main')
  if( IO_L ) write(IO_FID_LOG,*) '++++++ END   MAKING INITIAL DATA ++++++'
  if( IO_L ) write(IO_FID_LOG,*)

  !########## Finalize ##########

  call TIME_rapreport

  call FIO_finalize
  call MONIT_finalize
  ! stop MPI
  call PRC_MPIstop

  stop
  !=============================================================================
end program scaleinit
