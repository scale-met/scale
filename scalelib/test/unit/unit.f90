program unit
  use scale_precision
  use scale_stdio
  use scale_grid_index
  use scale_tracer

  use scale_process, only: &
     PRC_setup,    &
     PRC_MPIstart, &
     PRC_MPIsetup, &
     PRC_MPIfinish
  use scale_const, only: &
     CONST_setup
  use scale_comm, only: &
     COMM_setup
  use scale_grid, only: &
     DZ, DX, DY, &
     GRID_allocate, &
     GRID_generate

  use test_atmos_phy_tb_smg

  use test_atmos_dyn

  use test_comm

  character(len=H_MID), parameter :: MODELNAME = "Unit test"

  ! start MPI
  call PRC_MPIstart

  ! setup standard I/O
  call IO_setup( MODELNAME, .false. )

  ! setup MPI
  call PRC_MPIsetup( .false. )

  ! setup process
  call PRC_setup

  ! setup constants
  call CONST_setup

  call GRID_INDEX_setup

  TRACER_TYPE = 'SN14'
  call TRACER_setup

  ! setup horizontal/veritical grid system
  DZ = 500.0_RP
  DX = 500.0_RP
  DY = 500.0_RP
  call GRID_allocate
  call GRID_generate

  ! setup mpi communication
  call COMM_setup

  write(*,*) "test_comm"
  call test_comm_run

  write(*,*) "test_atmos_phy_tb_smg_run"
  call test_atmos_phy_tb_smg_run

  write(*,*) "test_atmos_dyn_run"
  call test_atmos_dyn_run

  call PRC_MPIfinish

end program unit
