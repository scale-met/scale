program unit
  use scale
  use scale_atmos_grid_cartesC_index
  use scale_tracer
  use scale_rm_process, only: &
     PRC_setup
  use scale_comm, only: &
     COMM_setup
  use scale_atmos_grid_cartesC, only: &
     DZ, DX, DY, &
     GRID_CARTESC_allocate, &
     GRID_CARTESC_generate
  use scale_atmos_hydrometeor, only: &
     ATMOS_HYDROMETEOR_regist

  use test_atmos_phy_tb_smg
  use test_atmos_dyn
  use test_comm
  implicit none

  character(len=H_MID), parameter :: APPNAME = "Unit test"
  integer :: q0

  ! scale setup
  call SCALE_init( APPNAME )

  ! setup process
  call PRC_setup

  call GRID_CARTESC_INDEX_setup

  call ATMOS_HYDROMETEOR_regist(q0, 1, 1, 0, (/'QV','QC'/), (/'QV','QC'/), (/"kg/kg","kg/kg"/) )

  ! setup horizontal/veritical grid system
  DZ = 500.0_RP
  DX = 500.0_RP
  DY = 500.0_RP
  call GRID_CARTESC_allocate
  call GRID_CARTESC_generate

  ! setup mpi communication
  call COMM_setup

  write(*,*) "test_comm"
  call test_comm_run

  write(*,*) "test_atmos_phy_tb_smg_run"
  call test_atmos_phy_tb_smg_run

  write(*,*) "test_atmos_dyn_run"
  call test_atmos_dyn_run

  call SCALE_finalize

end program unit
