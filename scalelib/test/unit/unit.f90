program unit
  use scale
  use scale_atmos_grid_cartesC_index
  use scale_tracer
  use scale_prc_cartesC, only: &
     PRC_CARTESC_setup, &
     PRC_CARTESC_finalize
  use scale_comm_cartesC, only: &
     COMM_setup, &
     COMM_regist, &
     COMM_finalize
  use scale_atmos_grid_cartesC, only: &
     ATMOS_GRID_CARTESC_allocate, &
     ATMOS_GRID_CARTESC_generate, &
     ATMOS_GRID_CARTESC_finalize
  use scale_atmos_hydrometeor, only: &
     ATMOS_HYDROMETEOR_setup, &
     ATMOS_HYDROMETEOR_regist, &
     ATMOS_HYDROMETEOR_finalize

  use test_atmos_phy_tb_smg
  use test_atmos_dyn
  use test_comm
  implicit none

  character(len=H_MID), parameter :: APPNAME = "Unit test"

  integer :: gid
  integer :: q0

  ! scale setup
  call SCALE_init( APPNAME )

  ! setup process
  call PRC_CARTESC_setup

  call ATMOS_GRID_CARTESC_INDEX_setup( KMAX=10, IMAX=10, JMAX=2, IBLOCK=5, JBLOCK=1 )

  call ATMOS_HYDROMETEOR_setup
  call ATMOS_HYDROMETEOR_regist( 1, 0, &
                                 (/'QV','QC'/), (/'QV','QC'/), (/"kg/kg","kg/kg"/), &
                                 q0 )

  ! setup mpi communication
  call COMM_setup
  call COMM_regist( KA, IA, JA, IHALO, JHALO, gid )

  ! setup horizontal/veritical grid system
  call ATMOS_GRID_CARTESC_allocate
  call ATMOS_GRID_CARTESC_generate( DZ=500.0_RP, DX=500.0_RP, DY=500.0_RP )

  write(*,*) "test_comm"
  call test_comm_run

  write(*,*) "test_atmos_phy_tb_smg_run"
  call test_atmos_phy_tb_smg_run

  write(*,*) "test_atmos_dyn_run"
  call test_atmos_dyn_run

  call ATMOS_GRID_CARTESC_finalize

  call COMM_finalize

  call ATMOS_HYDROMETEOR_finalize

  call PRC_CARTESC_finalize

  call SCALE_finalize

end program unit
