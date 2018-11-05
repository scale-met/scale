module test_atmos_phy_tb_smg

  !-----------------------------------------------------------------------------
  use scale_precision
  use scale_atmos_grid_cartesC_index
  use scale_tracer
  use dc_test, only: &
     AssertEqual, &
     AssertLessThan
  use scale_prc, only: &
     PRC_MPIbarrier
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  public :: test_atmos_phy_tb_smg_run

  !-----------------------------------------------------------------------------
  real(RP), allocatable :: qflx_sgs_momz(:,:,:,:)
  real(RP), allocatable :: qflx_sgs_momx(:,:,:,:)
  real(RP), allocatable :: qflx_sgs_momy(:,:,:,:)
  real(RP), allocatable :: qflx_sgs_rhot(:,:,:,:)
  real(RP), allocatable :: qflx_sgs_rhoq(:,:,:,:,:)

  real(RP), allocatable :: nu_C(:,:,:) ! eddy viscosity (center)
  real(RP), allocatable :: Pr  (:,:,:) ! Prantle number
  real(RP), allocatable :: Ri  (:,:,:) ! Richardson number
  real(RP), allocatable :: N2  (:,:,:) ! Brunt-Vaisala frequency

  real(RP), allocatable :: MOMZ(:,:,:)
  real(RP), allocatable :: MOMX(:,:,:)
  real(RP), allocatable :: MOMY(:,:,:)
  real(RP), allocatable :: RHOT(:,:,:)
  real(RP), allocatable :: DENS(:,:,:)
  real(RP), allocatable :: QTRC(:,:,:,:)
  real(RP), allocatable :: MOMZ_t(:,:,:)
  real(RP), allocatable :: MOMX_t(:,:,:)
  real(RP), allocatable :: MOMY_t(:,:,:)
  real(RP), allocatable :: RHOT_t(:,:,:)
  real(RP), allocatable :: RHOQ_t(:,:,:,:)

  real(RP), allocatable :: GSQRT(:,:,:,:)
  real(RP), allocatable :: J13G(:,:,:,:)
  real(RP), allocatable :: J23G(:,:,:,:)
  real(RP) :: J33G

  real(RP), allocatable :: MAPF(:,:,:,:)
  real(RP), allocatable :: FZ3D(:,:,:)

  real(DP) :: dt


  real(RP), allocatable :: ZERO(:,:,:,:)

  integer :: KME ! end of main region

  integer :: k, i, j, iq
  character(len=17) :: message
  character(len=4) :: rankname
  !-----------------------------------------------------------------------------
contains

  subroutine test_atmos_phy_tb_smg_run
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_io, only: &
     H_SHORT
  use scale_atmos_phy_tb_smg, only: &
     ATMOS_PHY_TB_SMG_setup
  use scale_atmos_grid_cartesC, only: &
     CDZ  => ATMOS_GRID_CARTESC_CDZ, &
     CDX  => ATMOS_GRID_CARTESC_CDX, &
     CDY  => ATMOS_GRID_CARTESC_CDX, &
     FZ   => ATMOS_GRID_CARTESC_FZ,  &
     CZ   => ATMOS_GRID_CARTESC_CZ,  &
     CBFZ => ATMOS_GRID_CARTESC_CBFZ
  use scale_prc, only: &
     PRC_myrank
  !-----------------------------------------------------------------------------
  implicit none
  !-----------------------------------------------------------------------------
  !
  !++ parameters & variables
  !
  !=============================================================================
  real(RP) :: CZ3D(  KA,IA,JA)
  integer :: i, j

  write(rankname,'(A,I2.2,A)') "(", PRC_myrank, ")"

  ! allocate
  allocate( qflx_sgs_momz(KA,IA,JA,3) )
  allocate( qflx_sgs_momx(KA,IA,JA,3) )
  allocate( qflx_sgs_momy(KA,IA,JA,3) )
  allocate( qflx_sgs_rhot(KA,IA,JA,3) )
  allocate( qflx_sgs_rhoq(KA,IA,JA,3,QA) )

  allocate( nu_C(KA,IA,JA) )
  allocate( Pr  (KA,IA,JA) )
  allocate( Ri  (KA,IA,JA) )
  allocate( N2  (KA,IA,JA) )

  allocate( MOMZ(KA,IA,JA) )
  allocate( MOMX(KA,IA,JA) )
  allocate( MOMY(KA,IA,JA) )
  allocate( RHOT(KA,IA,JA) )
  allocate( DENS(KA,IA,JA) )
  allocate( QTRC(KA,IA,JA,QA) )
  allocate( MOMZ_t(KA,IA,JA) )
  allocate( MOMX_t(KA,IA,JA) )
  allocate( MOMY_t(KA,IA,JA) )
  allocate( RHOT_t(KA,IA,JA) )
  allocate( RHOQ_t(KA,IA,JA,QA) )

  allocate( GSQRT(KA,IA,JA,7) )
  allocate( J13G(KA,IA,JA,7) )
  allocate( J23G(KA,IA,JA,7) )

  allocate( MAPF(IA,JA,2,4) )
  allocate( FZ3D(0:KA,IA,JA) )

  allocate( ZERO(KA,IA,JA,3) )


  do j = JS-1, JE+1
  do i = IS-1, IE+1
     FZ3D(:,i,j) = FZ(:)
     CZ3D(:,i,j) = CZ(:)
  end do
  end do

  MAPF(:,:,:,:) = 1.0_RP

  !########## Initial setup ##########
  call ATMOS_PHY_TB_smg_setup( &
       FZ3D, CZ3D, CDX, CDY, MAPF )

  ZERO(:,:,:,:) = 0.0_RP

  KME = KE
  do k = KS+1, KE
     if ( CBFZ(k) > 0.0_RP ) then
        KME = k - 1
        exit
     end if
  end do

  GSQRT(:,:,:,:) = 1.0_RP
  J13G(:,:,:,:) = 0.0_RP
  J23G(:,:,:,:) = 0.0_RP
  J33G          = 1.0_RP

  dt = 1.0_RP

  !########## test ##########

  ! introduced lower limiter for S2 in SMG, so test_zero and test_constant fail
  call test_zero

  call test_constant

  call test_big

  call test_double

end subroutine test_atmos_phy_tb_smg_run
!=============================================================================


subroutine test_zero
  use scale_atmos_phy_tb_smg, only: &
     ATMOS_PHY_TB_smg
  use scale_atmos_grid_cartesC, only: &
     FDZ  => ATMOS_GRID_CARTESC_FDZ,  &
     RCDZ => ATMOS_GRID_CARTESC_RCDZ, &
     RFDZ => ATMOS_GRID_CARTESC_RFDZ, &
     CDX  => ATMOS_GRID_CARTESC_CDX,  &
     FDX  => ATMOS_GRID_CARTESC_FDX,  &
     CDY  => ATMOS_GRID_CARTESC_CDY,  &
     FDY  => ATMOS_GRID_CARTESC_FDY
  call PRC_MPIbarrier
  write(*,*) rankname, "Test zero"

  MOMZ(:,:,:) = 0.0_RP
  MOMX(:,:,:) = 0.0_RP
  MOMY(:,:,:) = 0.0_RP
  RHOT(:,:,:) = 1.0_RP
  DENS(:,:,:) = 1.0_RP
  QTRC(:,:,:,:) = 0.0_RP

  N2(:,:,:) = 0.0_RP

  call ATMOS_PHY_TB_smg( &
       qflx_sgs_momz, qflx_sgs_momx, qflx_sgs_momy, & ! (out)
       qflx_sgs_rhot, qflx_sgs_rhoq,                & ! (out)
       MOMZ_t, MOMX_t, MOMY_t, RHOT_t, RHOQ_t,      & ! (out)
       nu_C, Ri, Pr,                                & ! (out)
       MOMZ, MOMX, MOMY, RHOT, DENS, QTRC, N2,      & ! (in)
       FZ3D, FDZ, RCDZ, RFDZ, CDX, FDX, CDY, FDY,   & ! (in)
       GSQRT, J13G, J23G, J33G, MAPF, dt            ) ! (in)

  call AssertEqual("qflx_sgs_momz", ZERO(KS:KE-1,IS:IE,JS:JE,:), qflx_sgs_momz(KS:KE-1,IS:IE,JS:JE,:), &
       significant_digits = RP*2-3, ignore_digits = -RP*4)
  call AssertEqual("qflx_sgs_momx", ZERO(KS:KE,IS:IE,JS:JE,:), qflx_sgs_momx(KS:KE,IS:IE,JS:JE,:), &
       significant_digits = RP*2-3, ignore_digits = -RP*4)
  call AssertEqual("qflx_sgs_momy", ZERO(KS:KE,IS:IE,JS:JE,:), qflx_sgs_momy(KS:KE,IS:IE,JS:JE,:), &
       significant_digits = RP*2-3, ignore_digits = -RP*4)
  call AssertEqual("qflx_sgs_rhot", ZERO(KS:KE,IS:IE,JS:JE,:), qflx_sgs_rhot(KS:KE,IS:IE,JS:JE,:), &
       significant_digits = RP*2-3, ignore_digits = -RP*4)
  call AssertEqual("MOMZ_t",        ZERO(KS:KE,IS:IE,JS:JE,1), MOMZ_t       (KS:KE,IS:IE,JS:JE), &
       significant_digits = RP*2-3, ignore_digits = -RP*3)
  call AssertEqual("MOMX_t",        ZERO(KS:KE,IS:IE,JS:JE,1), MOMX_t       (KS:KE,IS:IE,JS:JE), &
       significant_digits = RP*2-3, ignore_digits = -RP*3)
  call AssertEqual("MOMY_t",        ZERO(KS:KE,IS:IE,JS:JE,1), MOMY_t       (KS:KE,IS:IE,JS:JE), &
       significant_digits = RP*2-3, ignore_digits = -RP*3)
  call AssertEqual("RHOT_t",        ZERO(KS:KE,IS:IE,JS:JE,1), RHOT_t       (KS:KE,IS:IE,JS:JE), &
       significant_digits = RP*2-3, ignore_digits = -RP*3)
  do iq = 1, QA
     write(message, '("qflx_sgs_rhoq(",i2,")")') iq
     call AssertEqual(message, ZERO(KS:KE,IS:IE,JS:JE,:), qflx_sgs_rhoq(KS:KE,IS:IE,JS:JE,:,iq), &
       significant_digits = RP*2-3, ignore_digits = -RP*4)
     write(message, '("RHOQ_t(",i2,")")') iq
     call AssertEqual(message, ZERO(KS:KE,IS:IE,JS:JE,1), RHOQ_t       (KS:KE,IS:IE,JS:JE,  iq), &
       significant_digits = RP*2-3, ignore_digits = -RP*3)
  end do

end subroutine test_zero
!=============================================================================
subroutine test_constant
  use scale_atmos_phy_tb_smg, only: &
     ATMOS_PHY_TB_smg
  use scale_atmos_grid_cartesC, only: &
     FDZ  => ATMOS_GRID_CARTESC_FDZ,  &
     RCDZ => ATMOS_GRID_CARTESC_RCDZ, &
     RFDZ => ATMOS_GRID_CARTESC_RFDZ, &
     CDX  => ATMOS_GRID_CARTESC_CDX,  &
     FDX  => ATMOS_GRID_CARTESC_FDX,  &
     CDY  => ATMOS_GRID_CARTESC_CDY,  &
     FDY  => ATMOS_GRID_CARTESC_FDY

  call PRC_MPIbarrier
  write(*,*) rankname, "Test constant"

  MOMZ(:,:,:) = 1.0_RP
  MOMX(:,:,:) = 1.0_RP
  MOMY(:,:,:) = 1.0_RP
  RHOT(:,:,:) = 1.0_RP
  DENS(:,:,:) = 1.0_RP
  QTRC(:,:,:,:) = 1.0_RP

  N2(:,:,:) = 0.0_RP

  call fill_halo(MOMZ, MOMX, MOMY, RHOT, DENS, QTRC)

  call ATMOS_PHY_TB_smg( &
       qflx_sgs_momz, qflx_sgs_momx, qflx_sgs_momy, & ! (out)
       qflx_sgs_rhot, qflx_sgs_rhoq,                & ! (out)
       MOMZ_t, MOMX_t, MOMY_t, RHOT_t, RHOQ_t,      & ! (out)
       nu_C, Ri, Pr,                                & ! (out)
       MOMZ, MOMX, MOMY, RHOT, DENS, QTRC, N2,      & ! (in)
       FZ3D, FDZ, RCDZ, RFDZ, CDX, FDX, CDY, FDY,   & ! (in)
       GSQRT, J13G, J23G, J33G, MAPF, dt            ) ! (in)

  call AssertEqual("qflx_sgs_momz", ZERO(KS+1:KE-1,IS:IE,JS:JE,:), qflx_sgs_momz(KS+1:KE-1,IS:IE,JS:JE,:), &
       significant_digits = RP*2-3, ignore_digits = -RP*4)
  call AssertEqual("qflx_sgs_momx", ZERO(KS+1:KE-1,IS:IE,JS:JE,:), qflx_sgs_momx(KS+1:KE-1,IS:IE,JS:JE,:), &
       significant_digits = RP*2-3, ignore_digits = -RP*4)
  call AssertEqual("qflx_sgs_momy", ZERO(KS+1:KE-1,IS:IE,JS:JE,:), qflx_sgs_momy(KS+1:KE-1,IS:IE,JS:JE,:), &
       significant_digits = RP*2-3, ignore_digits = -RP*4)
  call AssertEqual("qflx_sgs_rhot", ZERO(KS:KE,IS:IE,JS:JE,:),     qflx_sgs_rhot(KS:KE,IS:IE,JS:JE,:), &
       significant_digits = RP*2-3, ignore_digits = -RP*4)
  call AssertEqual("MOMZ_t",        ZERO(KS+1:KE-2,IS:IE,JS:JE,1), MOMZ_t       (KS+1:KE-2,IS:IE,JS:JE), &
       significant_digits = RP*2-3, ignore_digits = -RP*3)
  call AssertEqual("MOMX_t",        ZERO(KS+2:KE-2,IS:IE,JS:JE,1), MOMX_t       (KS+2:KE-2,IS:IE,JS:JE), &
       significant_digits = RP*2-3, ignore_digits = -RP*3)
  call AssertEqual("MOMY_t",        ZERO(KS+2:KE-2,IS:IE,JS:JE,1), MOMY_t       (KS+2:KE-2,IS:IE,JS:JE), &
       significant_digits = RP*2-3, ignore_digits = -RP*3)
  call AssertEqual("RHOT_t",        ZERO(KS:KE,IS:IE,JS:JE,1),     RHOT_t       (KS:KE,IS:IE,JS:JE), &
       significant_digits = RP*2-3, ignore_digits = -RP*3)
  do iq = 1, QA
     write(message, '("qflx_sgs_rhoq(",i2,")")') iq
     call AssertEqual(message, ZERO(KS:KE,IS:IE,JS:JE,:), qflx_sgs_rhoq(KS:KE,IS:IE,JS:JE,:,iq), &
          significant_digits = RP*2-3, ignore_digits = -RP*4)
     write(message, '("RHOQ_t(",i2,")")') iq
     call AssertEqual(message, ZERO(KS:KE,IS:IE,JS:JE,1), RHOQ_t       (KS:KE,IS:IE,JS:JE,  iq), &
          significant_digits = RP*2-3, ignore_digits = -RP*3)
  end do

end subroutine test_constant
!=============================================================================
subroutine test_big
  use scale_const, only: &
     GRAV => CONST_GRAV
  use scale_atmos_grid_cartesC, only: &
     FDZ  => ATMOS_GRID_CARTESC_FDZ,  &
     RCDZ => ATMOS_GRID_CARTESC_RCDZ, &
     RFDZ => ATMOS_GRID_CARTESC_RFDZ, &
     CDX  => ATMOS_GRID_CARTESC_CDX,  &
     FDX  => ATMOS_GRID_CARTESC_FDX,  &
     CDY  => ATMOS_GRID_CARTESC_CDY,  &
     FDY  => ATMOS_GRID_CARTESC_FDY
  use scale_atmos_phy_tb_smg, only: &
     ATMOS_PHY_TB_smg

  real(RP) :: BIG(KA,IA,JA,3)

  real(RP) :: PI2

  PI2 = atan( 1.0_RP )*8.0_RP

  BIG(:,:,:,:) = 9.99E8_RP

  call PRC_MPIbarrier
  write(*,*) rankname, "Test big"
  ! check not to include BUG (UNDEF) value

  do j = 1, JA
  do i = 1, IA
  do k = 1, KA
     MOMZ(k,i,j) = 1.0_RP * sin( k*1.0_RP + i*2.0_RP + j*3.0_RP )
     MOMX(k,i,j) = 2.0_RP * cos( k*2.0_RP + i*3.0_RP + j*1.0_RP )
     MOMY(k,i,j) = 3.0_RP * sin( k*3.0_RP + i*1.0_RP + j*2.0_RP )
     RHOT(k,i,j) = 4.0_RP * cos( k*1.0_RP + i*1.0_RP + j*3.0_RP ) + 300.0_RP
     DENS(k,i,j) = 0.1_RP * sin( k*2.0_RP + i*2.0_RP + j*2.0_RP ) + 1.0_rp
     do iq = 1, QA
        QTRC(k,i,j,iq) = 6.0_RP * sin( k*3.0_RP + i*3.0_RP + j*1.0_RP + iq*2.0_RP ) + 6.0_RP
     end do
  end do
  end do
  end do

  do j = 1, JA
  do i = 1, IA
  do k = KS, KE
     N2(k,i,j) = GRAV * ( RHOT(k+1,i,j) - RHOT(k-1,i,j) ) / RHOT(k,i,j) * RCDZ(k)
  end do
  end do
  end do

  call fill_halo(MOMZ, MOMX, MOMY, RHOT, DENS, QTRC)

  call ATMOS_PHY_TB_smg( &
       qflx_sgs_momz, qflx_sgs_momx, qflx_sgs_momy, & ! (out)
       qflx_sgs_rhot, qflx_sgs_rhoq,                & ! (out)
       MOMZ_t, MOMX_t, MOMY_t, RHOT_t, RHOQ_t,      & ! (out)
       nu_C, Ri, Pr,                                & ! (out)
       MOMZ, MOMX, MOMY, RHOT, DENS, QTRC, N2,      & ! (in)
       FZ3D, FDZ, RCDZ, RFDZ, CDX, FDX, CDY, FDY,   & ! (in)
       GSQRT, J13G, J23G, J33G, MAPF, dt            ) ! (in)

  call AssertLessThan("qflx_sgs_momz", BIG(KS+1:KE-1,IS:IE,JS:JE,:), abs(qflx_sgs_momz(KS+1:KE-1,IS:IE,JS:JE,:)))
  call AssertLessThan("qflx_sgs_momx", BIG(KS:KE,IS:IE,JS:JE,:), abs(qflx_sgs_momx(KS:KE,IS:IE,JS:JE,:)))
  call AssertLessThan("qflx_sgs_momy", BIG(KS:KE,IS:IE,JS:JE,:), abs(qflx_sgs_momy(KS:KE,IS:IE,JS:JE,:)))
  call AssertLessThan("qflx_sgs_rhot", BIG(KS:KE,IS:IE,JS:JE,:), abs(qflx_sgs_rhot(KS:KE,IS:IE,JS:JE,:)))
  call AssertLessThan("MOMZ_t",        BIG(KS:KE,IS:IE,JS:JE,1), abs(MOMZ_t       (KS:KE,IS:IE,JS:JE  )))
  call AssertLessThan("MOMX_t",        BIG(KS:KE,IS:IE,JS:JE,1), abs(MOMX_t       (KS:KE,IS:IE,JS:JE  )))
  call AssertLessThan("MOMY_t",        BIG(KS:KE,IS:IE,JS:JE,1), abs(MOMY_t       (KS:KE,IS:IE,JS:JE  )))
  call AssertLessThan("RHOT_t",        BIG(KS:KE,IS:IE,JS:JE,1), abs(RHOT_t       (KS:KE,IS:IE,JS:JE  )))
  do iq = 1, QA
     if ( .not. TRACER_ADVC(iq) ) cycle
     write(message, '("qflx_sgs_rhoq(",i2,")")') iq
     call AssertLessThan(message, BIG(KS:KE,IS:IE,JS:JE,:), abs(qflx_sgs_rhoq(KS:KE,IS:IE,JS:JE,:,iq)))
     write(message, '("RHOQ_t(",i2,")")') iq
     call AssertLessThan(message, BIG(KS:KE,IS:IE,JS:JE,1), abs(RHOQ_t       (KS:KE,IS:IE,JS:JE,  iq)))
  end do

end subroutine test_big
!=============================================================================
subroutine test_double
  use scale_atmos_phy_tb_smg, only: &
     ATMOS_PHY_TB_smg
  use scale_atmos_grid_cartesC, only: &
     FDZ  => ATMOS_GRID_CARTESC_FDZ,  &
     RCDZ => ATMOS_GRID_CARTESC_RCDZ, &
     RFDZ => ATMOS_GRID_CARTESC_RFDZ, &
     CDX  => ATMOS_GRID_CARTESC_CDX,  &
     FDX  => ATMOS_GRID_CARTESC_FDX,  &
     CDY  => ATMOS_GRID_CARTESC_CDY,  &
     FDY  => ATMOS_GRID_CARTESC_FDY

  real(RP) :: qflx_sgs_momz2(KA,IA,JA,3)
  real(RP) :: qflx_sgs_momx2(KA,IA,JA,3)
  real(RP) :: qflx_sgs_momy2(KA,IA,JA,3)
  real(RP) :: qflx_sgs_rhot2(KA,IA,JA,3)
  real(RP) :: qflx_sgs_rhoq2(KA,IA,JA,3,QA)

  real(RP) :: work(KA,IA,JA,3)
  real(RP) :: FOUR(KA,IA,JA,3)
  real(RP) :: PI2

  FOUR(:,:,:,:) = 4.0_RP
  PI2 = atan( 1.0_RP )*8.0_RP

  call PRC_MPIbarrier
  write(*,*) rankname, "Test double"

  do j = 1, JA
  do i = 1, IA
  do k = 1, KA
     MOMZ(k,i,j) = 1.0_RP * sin( PI2 * ( k*1.0_RP/(KE-KS+1) + i*2.0_RP/(IE-IS+1) + j*3.0_RP/(JE-JS+1) ) )
     MOMX(k,i,j) = 2.0_RP * cos( PI2 * ( k*2.0_RP/(KE-KS+1) + i*3.0_RP/(IE-IS+1) + j*1.0_RP/(JE-JS+1) ) )
     MOMY(k,i,j) = 3.0_RP * sin( PI2 * ( k*3.0_RP/(KE-KS+1) + i*1.0_RP/(IE-IS+1) + j*2.0_RP/(JE-JS+1) ) )
     do iq = 1, QA
        QTRC(k,i,j,iq) = real(iq,RP) * sin( PI2 * ( k*2.0_RP/(KE-KS+1) + i*1.0_RP/(IE-IS+1) + j*3.0_RP/(JE-JS+1) ) )
     end do
  end do
  end do
  end do
  DENS(:,:,:) = 1.0_RP
  RHOT(:,:,:) = 1.0_RP ! Ri = 0
  N2(:,:,:) = 0.0_RP

  call fill_halo(MOMZ, MOMX, MOMY, RHOT, DENS, QTRC)

  call ATMOS_PHY_TB_smg( &
       qflx_sgs_momz, qflx_sgs_momx, qflx_sgs_momy, & ! (out)
       qflx_sgs_rhot, qflx_sgs_rhoq,                & ! (out)
       MOMZ_t, MOMX_t, MOMY_t, RHOT_t, RHOQ_t,      & ! (out)
       nu_C, Ri, Pr,                                & ! (out)
       MOMZ, MOMX, MOMY, RHOT, DENS, QTRC, N2,      & ! (in)
       FZ3D, FDZ, RCDZ, RFDZ, CDX, FDX, CDY, FDY,   & ! (in)
       GSQRT, J13G, J23G, J33G, MAPF, dt            ) ! (in)

  MOMZ(:,:,:) = MOMZ(:,:,:) * 2.0_RP
  MOMX(:,:,:) = MOMX(:,:,:) * 2.0_RP
  MOMY(:,:,:) = MOMY(:,:,:) * 2.0_RP
  QTRC(:,:,:,:) = QTRC(:,:,:,:) * 2.0_RP

  call ATMOS_PHY_TB_smg( &
       qflx_sgs_momz2, qflx_sgs_momx2, qflx_sgs_momy2, & ! (out)
       qflx_sgs_rhot2, qflx_sgs_rhoq2,                 & ! (out)
       MOMZ_t, MOMX_t, MOMY_t, RHOT_t, RHOQ_t,         & ! (out)
       nu_C, Ri, Pr,                                   & ! (out)
       MOMZ, MOMX, MOMY, RHOT, DENS, QTRC, N2,         & ! (in)
       FZ3D, FDZ, RCDZ, RFDZ, CDX, FDX, CDY, FDY,      & ! (in)
       GSQRT, J13G, J23G, J33G, MAPF, dt               ) ! (in)


  call AssertEqual("qflx_sgs_momz", FOUR(KS+1:KME-1,IS:IE,JS:JE,:), &
       qflx_sgs_momz2(KS+1:KME-1,IS:IE,JS:JE,:)/qflx_sgs_momz(KS+1:KME-1,IS:IE,JS:JE,:), &
       significant_digits = RP*2-3, ignore_digits = -RP*4)
  where(qflx_sgs_momx /= 0.0_RP) work = qflx_sgs_momx2 / qflx_sgs_momx
  where(qflx_sgs_momx == 0.0_RP) work = 4.0_RP
  call AssertEqual("qflx_sgs_momx", FOUR(KS:KE,IS:IE,JS:JE,:), work(KS:KE,IS:IE,JS:JE,:), &
       significant_digits = RP*2-3, ignore_digits = -RP*4)
  where(qflx_sgs_momy /= 0.0_RP) work = qflx_sgs_momy2 / qflx_sgs_momy
  where(qflx_sgs_momy == 0.0_RP) work = 4.0_RP
  call AssertEqual("qflx_sgs_momy", FOUR(KS:KE,IS:IE,JS:JE,:), work(KS:KE,IS:IE,JS:JE,:), &
       significant_digits = RP*2-3, ignore_digits = -RP*4)
  message = "iq = ??"
  do iq = 1, QA
     if ( .not. TRACER_ADVC(iq) ) cycle
     where(qflx_sgs_rhoq(:,:,:,:,iq) /= 0.0_RP) work = qflx_sgs_rhoq2(:,:,:,:,iq) / qflx_sgs_rhoq(:,:,:,:,iq)
     where(qflx_sgs_rhoq2(:,:,:,:,iq) == 0.0_RP) work = 4.0_RP
     write(message(6:7), "(i2)") iq
     call AssertEqual(message, FOUR(KS:KE,IS:IE,JS:JE,:), work(KS:KE,IS:IE,JS:JE,:), &
       significant_digits = RP*2-3, ignore_digits = -RP*4)
  end do

end subroutine test_double

!=============================================================================
subroutine fill_halo(MOMZ, MOMX, MOMY, RHOT, DENS, QTRC)
  implicit  none
  real(RP), intent(inout) :: MOMZ(KA,IA,JA)
  real(RP), intent(inout) :: MOMX(KA,IA,JA)
  real(RP), intent(inout) :: MOMY(KA,IA,JA)
  real(RP), intent(inout) :: RHOT(KA,IA,JA)
  real(RP), intent(inout) :: DENS(KA,IA,JA)
  real(RP), intent(inout) :: QTRC(KA,IA,JA,QA)

  do j = 1, JA
  do i = 1, IA
     MOMZ(1:KS-1,i,j) = MOMZ(KS,i,j)
     MOMX(1:KS-1,i,j) = MOMX(KS,i,j)
     MOMY(1:KS-1,i,j) = MOMY(KS,i,j)
     RHOT(1:KS-1,i,j) = RHOT(KS,i,j)
     DENS(1:KS-1,i,j) = DENS(KS,i,j)

     MOMZ(KE+1:KA,i,j) = MOMZ(KE,i,j)
     MOMX(KE+1:KA,i,j) = MOMX(KE,i,j)
     MOMY(KE+1:KA,i,j) = MOMY(KE,i,j)
     RHOT(KE+1:KA,i,j) = RHOT(KE,i,j)
     DENS(KE+1:KA,i,j) = DENS(KE,i,j)
  end do
  end do
  do iq = 1, QA
  do j = 1, JA
  do i = 1, IA
     QTRC(1:KS-1,i,j,iq) = QTRC(KS,i,j,iq)
     QTRC(KE+1:KA,i,j,iq) = QTRC(KE,i,j,iq)
  end do
  end do
  end do

end subroutine fill_halo


end module test_atmos_phy_tb_smg
