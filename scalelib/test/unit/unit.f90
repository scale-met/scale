program unit
  use mpi
  use scale_precision
  use scale_stdio
  use scale_grid_index
  use scale_tracer

  use scale_process, only: &
     PRC_UNIVERSAL_setup,  &
     PRC_global_setup,  &
     PRC_LOCAL_setup,  &
     PRC_MPIfinish, &
     PRC_MPIstart
  use scale_rm_process, only: &
     PRC_setup
  use scale_const, only: &
     CONST_setup
  use scale_comm, only: &
     COMM_setup
  use scale_grid, only: &
     DZ, DX, DY, &
     GRID_allocate, &
     GRID_generate
  use scale_atmos_hydrometeor, only: &
     ATMOS_HYDROMETEOR_regist

  use test_atmos_phy_tb_smg
  use test_atmos_dyn
  use test_comm
  implicit none

  character(len=H_MID), parameter :: MODELNAME = "Unit test"
  integer :: q0
  integer :: comm, myrank, nprocs
  logical :: ismaster

  ! start MPI
  call PRC_MPIstart( comm )

  ! setup standard I/O
  call IO_setup( MODELNAME, .false. )

  ! setup MPI
  call PRC_UNIVERSAL_setup( MPI_COMM_WORLD, nprocs, ismaster )
  call PRC_GLOBAL_setup( .true., MPI_COMM_WORLD )
  call PRC_LOCAL_setup( comm, myrank, ismaster )

  ! setup Log
  call IO_LOG_setup( myrank, ismaster )

  ! setup process
  call PRC_setup

  ! setup constants
  call CONST_setup

  call GRID_INDEX_setup

  call ATMOS_HYDROMETEOR_regist(q0, 1, 1, 0, (/'QV','QC'/), (/'QV','QC'/), (/"kg/kg","kg/kg"/) )

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
