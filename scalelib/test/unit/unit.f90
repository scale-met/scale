program unit
  use scale_stdio
  use scale_grid_index
  use scale_tracer

  use scale_process, only: &
     PRC_setup,    &
     PRC_MPIstart, &
     PRC_MPIfinish
  use scale_const, only: &
     CONST_setup
  use scale_comm, only: &
     COMM_setup
  use scale_grid, only: &
     GRID_allocate, &
     GRID_generate

  use test_atmos_phy_tb_smg

  use test_atmos_dyn_fent_fct

  character(len=H_MID), parameter :: MODELNAME = "Unit test"

  ! setup standard I/O
  call IO_setup( MODELNAME )

  ! start MPI
  call PRC_MPIstart

  ! setup process
  call PRC_setup

  ! setup constants
  call CONST_setup

  KMAX = 10
  IMAX = 10
  JMAX = 2
  IBLOCK = 5
  JBLOCK = 1

  call GRID_INDEX_setup

  TRACER_TYPE = 'SN14'
  call TRACER_setup

  ! setup horizontal/veritical grid system
  call GRID_allocate
  call GRID_generate

  ! setup mpi communication
  call COMM_setup

  write(*,*) "test_atmos_phy_tb_smg_run"
  call test_atmos_phy_tb_smg_run

  write(*,*) "test_atmos_dyn_fent_fct_run"
  call test_atmos_dyn_fent_fct_run

  call PRC_MPIfinish

end program unit
