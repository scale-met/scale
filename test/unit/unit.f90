program unit

  use mod_stdio, only: &
     IO_setup
  use mod_process, only: &
     PRC_setup,    &
     PRC_MPIstart, &
     PRC_MPIstop
  use mod_const, only: &
     CONST_setup
  use mod_comm, only: &
     COMM_setup
  use mod_grid, only: &
     GRID_allocate, &
     GRID_generate

  use test_atmos_phy_tb_smg

  use test_atmos_dyn_fent_fct

  ! setup standard I/O
  call IO_setup

  ! start MPI
  call PRC_MPIstart

  ! setup process
  call PRC_setup

  ! setup constants
  call CONST_setup

  ! setup horisontal/veritical grid system
  call GRID_allocate
  call GRID_generate

  ! setup mpi communication
  call COMM_setup

  call test_atmos_phy_tb_smg_run

  call test_atmos_dyn_fent_fct_run

  call PRC_MPIstop

end program unit
