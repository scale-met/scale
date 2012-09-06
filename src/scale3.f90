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
!! @li      2011-11-11 (H.Yashiro)  [new]
!! @li      2012-01-10 (H.Yashiro)  [mod] Change setup order, HISTORY module
!! @li      2012-03-23 (H.Yashiro)  [mod] GEOMETRICS, MONITOR module
!!
!<
!-------------------------------------------------------------------------------
#include "profiler.h"
  
program scaleles3
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_stdio, only: &
     IO_setup,   &
     IO_FID_CONF, &
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
     TIME_DOOCEAN_restart, &
     TIME_DOend,           &
     TIME_rapstart,        &
     TIME_rapend,          &
     TIME_rapreport,       &
     TIME_NSTEP_ATMOS_DYN     
  use mod_atmos_dyn_fent_fct_perf, only: &
       rk_ops, rk_min_ops
  use mod_perf, only: &
       rk_elapsed_time
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
  use mod_ocean_vars, only: &
     OCEAN_vars_restart_write, &
     OCEAN_sw_restart
  use dc_log, only: &
     LogInit
  use gtool_file, only: &
     FileCloseAll
  !-----------------------------------------------------------------------------
  implicit none
  !-----------------------------------------------------------------------------
  !
  !++ parameters & variables
  !
  !---------------------------------------------------------------------------
  integer*8 fp_ops, ld_ops, st_ops, ld_min_ops, st_min_ops
  integer num_iter
  !-----------------------------------------------------------------------------
  PROFILE_REGION_DECLARE(MAIN)
  PROFILE_REGION_DECLARE(INIT)
  !=============================================================================

  !########## Initial setup ##########

  ! setup standard I/O
  call IO_setup

  ! setup Log
  call LogInit(IO_FID_CONF, IO_FID_LOG, IO_L)

  ! start MPI
  call PRC_MPIstart

  ! setup process
  call PRC_setup

  ! setup constants
  call CONST_setup

  ! setup time
  call TIME_setup
  call TIME_rapstart('Initialize')
#ifdef _FPCOLL_
  call START_COLLECTION("Initialize")
#endif

  PROFILE_REGION_BEGIN(INIT)

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

  ! setup atmosphere
  call ATMOS_setup

  ! setup ocean
  call OCEAN_setup

#ifdef _FPCOLL_
  call STOP_COLLECTION  ("Initialize")
#endif
  call TIME_rapend('Initialize')

  PROFILE_REGION_END(INIT)

  num_iter = 0
  !########## main ##########

  if( IO_L ) write(IO_FID_LOG,*)
  if( IO_L ) write(IO_FID_LOG,*) '++++++ START TIMESTEP ++++++'
  call TIME_rapstart('Main Loop(Total)')
#ifdef _FPCOLL_
  call START_COLLECTION("Main")
#endif
  
  PROFILE_REGION_BEGIN(MAIN)
  
  do

    ! report current time
    call TIME_checkstate

    ! change to next state
    if ( TIME_DOATMOS_step ) call ATMOS_step
    if ( TIME_DOOCEAN_step ) call OCEAN_step

    ! time advance
    call TIME_advance

    ! history&monitor file output
    call HIST_write
    call MONIT_write('MAIN')

    ! restart output
    if ( ATMOS_sw_restart .AND. TIME_DOATMOS_restart ) call ATMOS_vars_restart_write
    if ( OCEAN_sw_restart .AND. TIME_DOOCEAN_restart ) call OCEAN_vars_restart_write

    num_iter = num_iter + 1
    if ( TIME_DOend ) exit

  enddo

#ifdef _FPCOLL_
  call STOP_COLLECTION  ("Main")
#endif

  PROFILE_REGION_END(MAIN)

  call TIME_rapend('Main Loop(Total)')
  if( IO_L ) write(IO_FID_LOG,*) '++++++ END TIMESTEP ++++++'
  if( IO_L ) write(IO_FID_LOG,*)


  call rk_ops(fp_ops, ld_ops, st_ops)
  call rk_min_ops(ld_min_ops, st_min_ops)
  fp_ops = fp_ops * TIME_NSTEP_ATMOS_DYN * num_iter
  ld_ops = ld_ops * TIME_NSTEP_ATMOS_DYN * num_iter
  st_ops = st_ops * TIME_NSTEP_ATMOS_DYN * num_iter
  ld_min_ops = ld_min_ops * TIME_NSTEP_ATMOS_DYN * num_iter
  st_min_ops = st_min_ops * TIME_NSTEP_ATMOS_DYN * num_iter
  write (*, '(A, F30.1, A)') "RK FLOPS: ", &
       fp_ops / rk_elapsed_time * 1.0e-9, " (GFLOPS)"
  write (*, '(A, F25.1, A)') "RK Throughput: ", &
       (ld_ops + st_ops) * 8 / rk_elapsed_time * 1.0e-9, " (GB/s)"
  write (*, '(A, F15.1, A)') "RK Effective Throughput: ", &
       (ld_min_ops + st_min_ops) * 8 / rk_elapsed_time * 1.0e-9, " (GB/s)"  

  !########## Finalize ##########

  call TIME_rapstart('Checkdiff')

  ! check data
  if ( ATMOS_sw_check ) call ATMOS_vars_restart_check

  call TIME_rapend('Checkdiff')

  call TIME_rapreport

  call FileCloseAll

  call MONIT_finalize
  ! stop MPI
  call PRC_MPIstop

  stop
  !=============================================================================
end program scaleles3
