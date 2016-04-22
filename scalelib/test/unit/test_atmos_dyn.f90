module test_atmos_dyn

  !-----------------------------------------------------------------------------
  use scale_precision
  use scale_grid_index
  use scale_index
  use scale_tracer
  use scale_atmos_dyn, only: &
     ATMOS_DYN
  use dc_test, only: &
     AssertEqual, &
     AssertGreaterThan, &
     AssertLessThan
  use scale_stdio, only: &
     H_SHORT
  use scale_comm, only: &
     COMM_vars8, &
     COMM_wait
  use scale_grid, only: &
       CZ   => GRID_CZ,   &
       FZ   => GRID_FZ,   &
       CDZ  => GRID_CDZ,  &
       CDX  => GRID_CDX,  &
       CDY  => GRID_CDY,  &
       FDZ  => GRID_FDZ,  &
       FDX  => GRID_FDX,  &
       FDY  => GRID_FDY,  &
       RCDZ => GRID_RCDZ, &
       RCDX => GRID_RCDX, &
       RCDY => GRID_RCDY, &
       RFDZ => GRID_RFDZ, &
       RFDX => GRID_RFDX, &
       RFDY => GRID_RFDY
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  public :: test_atmos_dyn_run

  !-----------------------------------------------------------------------------
  real(RP), allocatable :: DENS(:,:,:)
  real(RP), allocatable :: MOMZ(:,:,:)
  real(RP), allocatable :: MOMX(:,:,:)
  real(RP), allocatable :: MOMY(:,:,:)
  real(RP), allocatable :: RHOT(:,:,:)
  real(RP), allocatable :: QTRC(:,:,:,:)

  real(RP), allocatable :: DENS_av(:,:,:)
  real(RP), allocatable :: MOMZ_av(:,:,:)
  real(RP), allocatable :: MOMX_av(:,:,:)
  real(RP), allocatable :: MOMY_av(:,:,:)
  real(RP), allocatable :: RHOT_av(:,:,:)
  real(RP), allocatable :: QTRC_av(:,:,:,:)

  real(RP), allocatable :: DENS_tp(:,:,:)
  real(RP), allocatable :: MOMZ_tp(:,:,:)
  real(RP), allocatable :: MOMX_tp(:,:,:)
  real(RP), allocatable :: MOMY_tp(:,:,:)
  real(RP), allocatable :: RHOT_tp(:,:,:)
  real(RP), allocatable :: QTRC_tp(:,:,:,:)

  real(RP), allocatable :: DENS_o(:,:,:)
  real(RP), allocatable :: MOMZ_o(:,:,:)
  real(RP), allocatable :: MOMX_o(:,:,:)
  real(RP), allocatable :: MOMY_o(:,:,:)
  real(RP), allocatable :: RHOT_o(:,:,:)
  real(RP), allocatable :: QTRC_o(:,:,:,:)

  real(RP), allocatable :: REF_dens(:,:,:)
  real(RP), allocatable :: REF_pott(:,:,:)
  real(RP), allocatable :: REF_qv  (:,:,:)
  real(RP), allocatable :: REF_pres(:,:,:)

  real(RP), allocatable :: DAMP_var(:,:,:,:)
  real(RP), allocatable :: DAMP_alpha(:,:,:,:)

  real(RP), allocatable :: PHI(:,:,:)
  real(RP), allocatable :: GSQRT(:,:,:,:)
  real(RP), allocatable :: J13G(:,:,:,:)
  real(RP), allocatable :: J23G(:,:,:,:)
  real(RP) :: J33G
  real(RP), allocatable :: MAPF(:,:,:,:)

  real(RP), allocatable :: AQ_CV(:)

  integer  :: nd_order
  real(RP) :: nd_coef
  real(RP) :: nd_sfc_fact
  logical  :: nd_use_rs

  real(RP) :: divdmp_coef

  logical  :: flag_fct_rho      = .true.
  logical  :: flag_fct_momentum = .true.
  logical  :: flag_fct_t        = .true.
  logical  :: flag_fct_along_stream = .true.

  real(RP), allocatable :: ZERO(:,:,:)

  integer :: KME ! end of main region

  character(len=H_SHORT) :: DYN_TYPE

  integer :: k, i, j, iq
  character(len=11) :: message
  !-----------------------------------------------------------------------------
contains

  subroutine test_atmos_dyn_run
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_atmos_dyn, only: &
     ATMOS_DYN_setup
  use scale_grid, only: &
     GRID_CBFZ
  use scale_const, only: &
     GRAV => CONST_GRAV

  !-----------------------------------------------------------------------------
  implicit none
  !-----------------------------------------------------------------------------
  !
  !++ parameters & variables
  !
  !-----------------------------------------------------------------------------
  real(RP) :: lat(1,IA,JA)
  integer :: j
  !=============================================================================

  !########## Initial setup ##########

  ! allocate
  allocate( DENS(KA,IA,JA) )
  allocate( MOMZ(KA,IA,JA) )
  allocate( MOMX(KA,IA,JA) )
  allocate( MOMY(KA,IA,JA) )
  allocate( RHOT(KA,IA,JA) )
  allocate( QTRC(KA,IA,JA,QA) )

  allocate( DENS_av(KA,IA,JA) )
  allocate( MOMZ_av(KA,IA,JA) )
  allocate( MOMX_av(KA,IA,JA) )
  allocate( MOMY_av(KA,IA,JA) )
  allocate( RHOT_av(KA,IA,JA) )
  allocate( QTRC_av(KA,IA,JA,QA) )

  allocate( DENS_tp(KA,IA,JA) )
  allocate( MOMZ_tp(KA,IA,JA) )
  allocate( MOMX_tp(KA,IA,JA) )
  allocate( MOMY_tp(KA,IA,JA) )
  allocate( RHOT_tp(KA,IA,JA) )
  allocate( QTRC_tp(KA,IA,JA,QA) )

  allocate( DENS_o(KA,IA,JA) )
  allocate( MOMZ_o(KA,IA,JA) )
  allocate( MOMX_o(KA,IA,JA) )
  allocate( MOMY_o(KA,IA,JA) )
  allocate( RHOT_o(KA,IA,JA) )
  allocate( QTRC_o(KA,IA,JA,QA) )

  allocate( REF_dens(KA,IA,JA) )
  allocate( REF_pott(KA,IA,JA) )
  allocate( REF_qv  (KA,IA,JA) )
  allocate( REF_pres(KA,IA,JA) )

  allocate( DAMP_var(KA,IA,JA,5+QA) )
  allocate( DAMP_alpha(KA,IA,JA,5+QA) )

  allocate( PHI(KA,IA,JA) )
  allocate( GSQRT(KA,IA,JA,7) )
  allocate( J13G(KA,IA,JA,4) )
  allocate( J23G(KA,IA,JA,4) )

  allocate( MAPF(IA,JA,2,4) )

  allocate( AQ_CV(QA) )

  allocate( ZERO(KA,IA,JA) )

  ZERO(:,:,:) = 0.0_RP

  nd_order = 2
  nd_coef = 0.01_RP
  nd_sfc_fact = 1.0_RP
  nd_use_rs = .true.
  do j = 1, JA
     lat(1,:,j) = real(j, RP)
  end do

  DYN_TYPE = "HEVE"
  call ATMOS_DYN_setup( DYN_TYPE,                           & ! (in)
                        DENS, MOMZ, MOMX, MOMY, RHOT, QTRC, & ! (in)
                        CDZ, CDX, CDY, FDZ, FDX, FDY,       & ! (in)
                        .false., lat                        ) ! (in)

  do k = KS+1, KE
     if ( GRID_CBFZ(k) > 0.0_RP ) then
        KME = k - 1
        exit
     end if
  end do

  MOMZ(KE,:,:) = 0.0_RP

  AQ_CV(:) = 1.0_RP

  do j = 1, JA
  do i = 1, IA
     PHI(:,i,j) = CZ(:) * GRAV
  end do
  end do
  GSQRT(:,:,:,:) = 1.0_RP
  J13G(:,:,:,:) = 0.0_RP
  J23G(:,:,:,:) = 0.0_RP
  J33G          = 1.0_RP
  MAPF(:,:,:,:) = 1.0_RP

  divdmp_coef = 0.0_RP

  DENS_tp(:,:,:) = 0.0_RP
  MOMZ_tp(:,:,:) = 0.0_RP
  MOMX_tp(:,:,:) = 0.0_RP
  MOMY_tp(:,:,:) = 0.0_RP
  RHOT_tp(:,:,:) = 0.0_RP
  QTRC_tp(:,:,:,:) = 0.0_RP

  !########## test ##########


  call test_undef

  call test_const

  call test_conserve

  call test_cwc

  call test_fctminmax

end subroutine test_atmos_dyn_run
!=============================================================================


subroutine test_undef

  real(RP) :: BIG(KA,IA,JA)

  BIG(:,:,:) = 9.99E9_RP

  write(*,*) "Test undef"

  do j = 1, JA
  do i = 1, IA
  do k = 1, KA
     MOMZ(k,i,j) = k + i + j
     MOMX(k,i,j) = k + i - j
     MOMY(k,i,j) = k - i + j
     RHOT(k,i,j) = k - i - j + KA + IA + JA
     DENS(k,i,j) = k + i + j + 1.0_RP
     do iq = 1, QA
        QTRC(k,i,j,iq) = 1.0_RP*(k + i + j + iq) / (KA + IA + JA + QA) / QA
     end do
  end do
  end do
  end do

  do j = 1, JA
  do i = 1, IA
  do k = 1, KA
     REF_dens(k,i,j) = KA - k + 1
     REF_pott(k,i,j) = k
     REF_qv  (k,i,j) = KA - k + 1
     REF_pres(k,i,j) = KA - k + 1
  end do
  end do
  end do

  DAMP_var  (:,:,:,:) = -9.999E30_RP
  DAMP_alpha(:,:,:,:) = 0.0_RP

  do i = 1, 2
     call ATMOS_DYN( &
          DENS, MOMZ, MOMX, MOMY, RHOT, QTRC,          & ! (inout)
          DENS_av, MOMZ_av, MOMX_av, MOMY_av, RHOT_av, QTRC_av, & ! (out)
          DENS_tp, MOMZ_tp, MOMX_tp, MOMY_tp, RHOT_tp, QTRC_tp, & ! (in)
          CDZ, CDX, CDY, FDZ, FDX, FDY,                & ! (in)
          RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,          & ! (in)
          PHI, GSQRT, J13G, J23G, J33G, MAPF,          & ! (in)
          AQ_CV,                                       & ! (in)
          REF_dens, REF_pott, REF_qv, REF_pres,        & ! (in)
          nd_coef, nd_coef, nd_order, nd_sfc_fact, nd_use_rs, & ! (in)
          DAMP_var(:,:,:,1), DAMP_var(:,:,:,2), DAMP_var(:,:,:,3), DAMP_var(:,:,:,4), DAMP_var(:,:,:,5), DAMP_var(:,:,:,6:6+QA-1), & ! (in)
          DAMP_alpha(:,:,:,1), DAMP_alpha(:,:,:,2), DAMP_alpha(:,:,:,3), DAMP_alpha(:,:,:,4), DAMP_alpha(:,:,:,5), & ! (in)
          DAMP_alpha(:,:,:,6:6+QA-1),                  & ! (in)
          divdmp_coef,                                 & ! (in)
          flag_fct_rho, flag_fct_momentum, flag_fct_t, & ! (in)
          flag_fct_along_stream,                       & ! (in)
          .false.,                                     & ! (in)
          1.0_DP, 1.0_DP, 1                            ) ! (in)
  end do

  call AssertLessThan("MOMZ", BIG(KS:KE,IS:IE,JS:JE), MOMZ(KS:KE,IS:IE,JS:JE))
  call AssertLessThan("MOMX", BIG(KS:KE,IS:IE,JS:JE), MOMX(KS:KE,IS:IE,JS:JE))
  call AssertLessThan("MOMY", BIG(KS:KE,IS:IE,JS:JE), MOMY(KS:KE,IS:IE,JS:JE))
  call AssertLessThan("RHOT", BIG(KS:KE,IS:IE,JS:JE), RHOT(KS:KE,IS:IE,JS:JE))
  message = "iq = ??"
  do iq = 1, QA
     write(message(6:7), "(i2)") iq
     call AssertLessThan(message, BIG(KS:KE,IS:IE,JS:JE), QTRC(KS:KE,IS:IE,JS:JE,iq))
  end do

end subroutine test_undef

subroutine test_const
  real(RP) :: answer(KA,IA,JA)

  write(*,*) "Test constant"

  DENS(:,:,:) = 1.0_RP
  MOMZ(:,:,:) = 2.0_RP
  MOMX(:,:,:) = 3.0_RP
  MOMY(:,:,:) = 4.0_RP
  RHOT(:,:,:) = 300.0_RP
  QTRC(:,:,:,:) = 0.1_RP

  REF_dens(:,:,:) = 1.0_RP
  REF_pott(:,:,:) = 300.0_RP
  REF_qv(:,:,:)   = 0.001_RP
  REF_pres(:,:,:) = 1000._RP

  DAMP_var  (:,:,:,:) = -9.999E30_RP
  DAMP_alpha(:,:,:,:) = 0.0_RP

  call ATMOS_DYN( &
       DENS, MOMZ, MOMX, MOMY, RHOT, QTRC,          & ! (inout)
       DENS_av, MOMZ_av, MOMX_av, MOMY_av, RHOT_av, QTRC_av, & ! (out)
       DENS_tp, MOMZ_tp, MOMX_tp, MOMY_tp, RHOT_tp, QTRC_tp, & ! (in)
       CDZ, CDX, CDY, FDZ, FDX, FDY,                & ! (in)
       RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,          & ! (in)
       PHI, GSQRT, J13G, J23G, J33G, MAPF,          & ! (in)
       AQ_CV,                                       & ! (in)
       REF_dens, REF_pott, REF_qv, REF_pres,        & ! (in)
       nd_coef, nd_coef, nd_order, nd_sfc_fact, nd_use_rs, & ! (in)
       DAMP_var(:,:,:,1), DAMP_var(:,:,:,2), DAMP_var(:,:,:,3), DAMP_var(:,:,:,4), DAMP_var(:,:,:,5), DAMP_var(:,:,:,6:6+QA-1), & ! (in)
       DAMP_alpha(:,:,:,1), DAMP_alpha(:,:,:,2), DAMP_alpha(:,:,:,3), DAMP_alpha(:,:,:,4), DAMP_alpha(:,:,:,5), & ! (in)
       DAMP_alpha(:,:,:,6:6+QA-1),                  & ! (in)
       divdmp_coef,                                 & ! (in)
       flag_fct_rho, flag_fct_momentum, flag_fct_t, & ! (in)
       flag_fct_along_stream,                       & ! (in)
       .false.,                                     & ! (in)
       1.0_DP, 1.0_DP, 1                            ) ! (in)

  do k = KS, KE
     answer(k,:,:) = MOMZ(k,IS,JS)
  end do
  call AssertEqual("MOMZ", answer(KS:KE,IS:IE,JS:JE), MOMZ(KS:KE,IS:IE,JS:JE))
  do k = KS, KE
     answer(k,:,:) = MOMX(k,IS,JS)
  end do
  call AssertEqual("MOMX", answer(KS:KE,IS:IE,JS:JE), MOMX(KS:KE,IS:IE,JS:JE))
  do k = KS, KE
     answer(k,:,:) = MOMY(k,IS,JS)
  end do
  call AssertEqual("MOMY", answer(KS:KE,IS:IE,JS:JE), MOMY(KS:KE,IS:IE,JS:JE))
  do k = KS, KE
     answer(k,:,:) = DENS(k,IS,JS)
  end do
  call AssertEqual("DENS", answer(KS:KE,IS:IE,JS:JE), DENS(KS:KE,IS:IE,JS:JE))
  do k = KS, KE
     answer(k,:,:) = RHOT(k,IS,JS)
  end do
  call AssertEqual("RHOT", answer(KS:KE,IS:IE,JS:JE), RHOT(KS:KE,IS:IE,JS:JE))
  message = "iq = ??"
  do iq = 1, QA
     do k = KS, KE
        answer(k,:,:) = QTRC(k,IS,JS,iq)
     end do
     write(message(6:7), "(i2)") iq
     call AssertEqual(message, answer(KS:KE,IS:IE,JS:JE), QTRC(KS:KE,IS:IE,JS:JE,iq))
  end do

end subroutine test_const

subroutine test_conserve

  real(RP) :: total_o, total

  write(*,*) "Test conserve"

  do j = 1, JA
  do i = 1, IA
  do k = KS, KE
     MOMZ(k,i,j) = (k + i + j) * 1.0_RP
     MOMX(k,i,j) = (k*2 + i + j) * 1.0_RP
     MOMY(k,i,j) = (k + i*2 + j) * 1.0_RP
     RHOT(k,i,j) = (k + i + j*2) * 1.0_RP
     DENS(k,i,j) = (k*2 + i*2 + j) * 1.0_RP + 1.0_RP
     do iq = 1, QA
        QTRC(k,i,j,iq) = (k + i + j + iq) * 0.1_RP / (KE + IE + JE + QA)
     end do
  end do
  end do
  end do

  REF_dens(:,:,:) = 1.0_RP
  REF_pott(:,:,:) = 1.0_RP
  REF_qv(:,:,:)   = 1.0_RP
  REF_pres(:,:,:) = 1.0_RP

  DAMP_var  (:,:,:,:) = -9.999E30_RP
  DAMP_alpha(:,:,:,:) = 0.0_RP

  call COMM_vars8( DENS(:,:,:), 1 )
  call COMM_vars8( MOMZ(:,:,:), 2 )
  call COMM_vars8( MOMX(:,:,:), 3 )
  call COMM_vars8( MOMY(:,:,:), 4 )
  call COMM_vars8( RHOT(:,:,:), 5 )
  call COMM_wait ( DENS(:,:,:), 1 )
  call COMM_wait ( MOMZ(:,:,:), 2 )
  call COMM_wait ( MOMX(:,:,:), 3 )
  call COMM_wait ( MOMY(:,:,:), 4 )
  call COMM_wait ( RHOT(:,:,:), 5 )

  do iq = 1, QA
     call COMM_vars8( QTRC(:,:,:,iq), iq )
  enddo
  do iq = 1, QA
     call COMM_wait ( QTRC(:,:,:,iq), iq )
  enddo

  call copy

  call ATMOS_DYN( &
         DENS, MOMZ, MOMX, MOMY, RHOT, QTRC,          & ! (out)
         DENS_av, MOMZ_av, MOMX_av, MOMY_av, RHOT_av, QTRC_av, & ! (inout)
         DENS_tp, MOMZ_tp, MOMX_tp, MOMY_tp, RHOT_tp, QTRC_tp, & ! (in)
         CDZ, CDX, CDY, FDZ, FDX, FDY,                & ! (in)
         RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,          & ! (in)
         PHI, GSQRT, J13G, J23G, J33G, MAPF,          & ! (in)
         AQ_CV,                                       & ! (in)
         REF_dens, REF_pott, REF_qv, REF_pres,        & ! (in)
         nd_coef, nd_coef, nd_order, nd_sfc_fact, nd_use_rs, & ! (in)
         DAMP_var(:,:,:,1), DAMP_var(:,:,:,2), DAMP_var(:,:,:,3), DAMP_var(:,:,:,4), DAMP_var(:,:,:,5), DAMP_var(:,:,:,6:6+QA-1), & ! (in)
         DAMP_alpha(:,:,:,1), DAMP_alpha(:,:,:,2), DAMP_alpha(:,:,:,3), DAMP_alpha(:,:,:,4), DAMP_alpha(:,:,:,5), & ! (in)
         DAMP_alpha(:,:,:,6:6+QA-1),                  & ! (in)
         divdmp_coef,                                 & ! (in)
         flag_fct_rho, flag_fct_momentum, flag_fct_t, & ! (in)
         flag_fct_along_stream,                       & ! (in)
         .true.,                                      & ! (in)
         1.0_RP, 1.0_DP, 1                            ) ! (in)

  total_o = 0.0_RP
  total = 0.0_RP
  do j = JS, JE
  do i = IS, IE
  do k = KS, KE
     total_o = total_o + DENS_o(k,i,j) * CDZ(k) * CDX(i) * CDY(j)
     total   = total   + DENS  (k,i,j) * CDZ(k) * CDX(i) * CDY(j)
  end do
  end do
  end do
  call AssertEqual("DENS", total_o, total, RP*2-2, -10)

  total_o = 0.0_RP
  total = 0.0_RP
  do j = JS, JE
  do i = IS, IE
  do k = KS, KE
     total_o = total_o + MOMX_o(k,i,j) * CDZ(k) * CDX(i) * CDY(j)
     total   = total   + MOMX  (k,i,j) * CDZ(k) * CDX(i) * CDY(j)
  end do
  end do
  end do
  call AssertEqual("MOMX", total_o, total, RP*2-2, -10)

  total_o = 0.0_RP
  total = 0.0_RP
  do j = JS, JE
  do i = IS, IE
  do k = KS, KE
     total_o = total_o + MOMY_o(k,i,j) * CDZ(k) * CDX(i) * CDY(j)
     total   = total   + MOMY  (k,i,j) * CDZ(k) * CDX(i) * CDY(j)
  end do
  end do
  end do
  call AssertEqual("MOMY", total_o, total, RP*2-2, -10)

  total_o = 0.0_RP
  total = 0.0_RP
  do j = JS, JE
  do i = IS, IE
  do k = KS, KE
     total_o = total_o + RHOT_o(k,i,j) * CDZ(k) * CDX(i) * CDY(j)
     total   = total   + RHOT  (k,i,j) * CDZ(k) * CDX(i) * CDY(j)
  end do
  end do
  end do
  call AssertEqual("RHOT", total_o, total, RP*2-2, -10)

  message = "iq = ??"
  do iq = 1, QA
     write(message(6:7), "(i2)") iq
     total_o = 0.0_RP
     total = 0.0_RP
     do j = JS, JE
     do i = IS, IE
     do k = KS, KE
        total_o = total_o + RHOT_o(k,i,j) * CDZ(k) * CDX(i) * CDY(j)
        total   = total   + RHOT  (k,i,j) * CDZ(k) * CDX(i) * CDY(j)
     end do
     end do
     end do
     call AssertEqual(message, total_o, total, RP*2-2, -10)
  end do

end subroutine test_conserve

subroutine test_cwc
  real(RP), parameter :: POTT = 300._RP
  real(RP), parameter :: Q = 0.01_RP
  real(RP) :: answer(KA,IA,JA)

  write(*,*) "Test CWC"

  call random_number(DENS)
  DENS(:,:,:) = DENS(:,:,:)*0.1 + 1.0_RP
  call random_number(MOMZ)
  call random_number(MOMX)
  call random_number(MOMY)

  RHOT(:,:,:) = POTT * DENS(:,:,:)
  QTRC(:,:,:,:) = Q

  REF_dens(:,:,:) = 1.0_RP
  REF_pott(:,:,:) = POTT
  REF_qv(:,:,:)   = Q
  REF_pres(:,:,:) = 1000._RP

  DAMP_var  (:,:,:,:) = -9.999E30_RP
  DAMP_alpha(:,:,:,:) = 0.0_RP

  call COMM_vars8( DENS(:,:,:), 1 )
  call COMM_vars8( MOMZ(:,:,:), 2 )
  call COMM_vars8( MOMX(:,:,:), 3 )
  call COMM_vars8( MOMY(:,:,:), 4 )
  call COMM_vars8( RHOT(:,:,:), 5 )
  call COMM_wait ( DENS(:,:,:), 1 )
  call COMM_wait ( MOMZ(:,:,:), 2 )
  call COMM_wait ( MOMX(:,:,:), 3 )
  call COMM_wait ( MOMY(:,:,:), 4 )
  call COMM_wait ( RHOT(:,:,:), 5 )

  do iq = 1, QA
     call COMM_vars8( QTRC(:,:,:,iq), iq )
  enddo
  do iq = 1, QA
     call COMM_wait ( QTRC(:,:,:,iq), iq )
  enddo

  call ATMOS_DYN( &
       DENS, MOMZ, MOMX, MOMY, RHOT, QTRC,          & ! (inout)
       DENS_av, MOMZ_av, MOMX_av, MOMY_av, RHOT_av, QTRC_av, & ! (out)
       DENS_tp, MOMZ_tp, MOMX_tp, MOMY_tp, RHOT_tp, QTRC_tp, & ! (in)
       CDZ, CDX, CDY, FDZ, FDX, FDY,                & ! (in)
       RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,          & ! (in)
       PHI, GSQRT, J13G, J23G, J33G, MAPF,          & ! (in)
       AQ_CV,                                       & ! (in)
       REF_dens, REF_pott, REF_qv, REF_pres,        & ! (in)
       nd_coef, nd_coef, nd_order, nd_sfc_fact, nd_use_rs, & ! (in)
       DAMP_var(:,:,:,1), DAMP_var(:,:,:,2), DAMP_var(:,:,:,3), DAMP_var(:,:,:,4), DAMP_var(:,:,:,5), DAMP_var(:,:,:,6:6+QA-1), & ! (in)
       DAMP_alpha(:,:,:,1), DAMP_alpha(:,:,:,2), DAMP_alpha(:,:,:,3), DAMP_alpha(:,:,:,4), DAMP_alpha(:,:,:,5), & ! (in)
       DAMP_alpha(:,:,:,6:6+QA-1),                  & ! (in)
       divdmp_coef,                                 & ! (in)
       flag_fct_rho, flag_fct_momentum, flag_fct_t, & ! (in)
       flag_fct_along_stream,                       & ! (in)
       .false.,                                     & ! (in)
       1.0_DP, 1.0_DP, 1                            ) ! (in)

  answer(:,:,:) = POTT
  call AssertEqual("POTT", answer(KS:KE,IS:IE,JS:JE), RHOT(KS:KE,IS:IE,JS:JE)/DENS(KS:KE,IS:IE,JS:JE), RP*2-2, -10)

  answer(:,:,:) = Q
  message = "iq = ??"
  do iq = 1, QA
     write(message(6:7), "(i2)") iq
     call AssertEqual(message, answer(KS:KE,IS:IE,JS:JE), QTRC(KS:KE,IS:IE,JS:JE,iq), RP*2-2, -10)
  end do

end subroutine test_cwc

subroutine test_fctminmax
  real(RP), parameter :: POTT_MAX = 350_RP
  real(RP), parameter :: POTT_MIN = 300_RP
  real(RP), parameter :: Q_MAX = 0.1_RP
  real(RP), parameter :: Q_MIN = 0_RP
  real(RP) :: MINMAX(KA,IA,JA)
  real(RP) :: epsilon
  integer :: k, i, j

  epsilon = 0.1_RP**(RP*2-2)

  write(*,*) "Test FCT"

  call random_number(DENS)
  DENS(:,:,:) = DENS(:,:,:)*0.1 + 1.0_RP
  call random_number(MOMZ)
  call random_number(MOMX)
  call random_number(MOMY)

  RHOT(:,:,:) = POTT_MIN * DENS(:,:,:)
  QTRC(:,:,:,:) = Q_MIN
  do j = JS, JS+JMAX/2
  do i = IS, IS+IMAX/2
  do k = KS, KS+KMAX/2
     RHOT(k,i,j) = POTT_MAX * DENS(k,i,j)
     QTRC(k,i,j,:) = Q_MAX
  end do
  end do
  end do

  REF_dens(:,:,:) = 1.0_RP
  REF_pott(:,:,:) = POTT_MIN
  REF_qv(:,:,:)   = Q_MIN
  REF_pres(:,:,:) = 1000._RP

  DAMP_var  (:,:,:,:) = -9.999E30_RP
  DAMP_alpha(:,:,:,:) = 0.0_RP

  call COMM_vars8( DENS(:,:,:), 1 )
  call COMM_vars8( MOMZ(:,:,:), 2 )
  call COMM_vars8( MOMX(:,:,:), 3 )
  call COMM_vars8( MOMY(:,:,:), 4 )
  call COMM_vars8( RHOT(:,:,:), 5 )
  call COMM_wait ( DENS(:,:,:), 1 )
  call COMM_wait ( MOMZ(:,:,:), 2 )
  call COMM_wait ( MOMX(:,:,:), 3 )
  call COMM_wait ( MOMY(:,:,:), 4 )
  call COMM_wait ( RHOT(:,:,:), 5 )

  do iq = 1, QA
     call COMM_vars8( QTRC(:,:,:,iq), iq )
  enddo
  do iq = 1, QA
     call COMM_wait ( QTRC(:,:,:,iq), iq )
  enddo

  call ATMOS_DYN( &
       DENS, MOMZ, MOMX, MOMY, RHOT, QTRC,          & ! (inout)
       DENS_av, MOMZ_av, MOMX_av, MOMY_av, RHOT_av, QTRC_av, & ! (out)
       DENS_tp, MOMZ_tp, MOMX_tp, MOMY_tp, RHOT_tp, QTRC_tp, & ! (in)
       CDZ, CDX, CDY, FDZ, FDX, FDY,                & ! (in)
       RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,          & ! (in)
       PHI, GSQRT, J13G, J23G, J33G, MAPF,          & ! (in)
       AQ_CV,                                       & ! (in)
       REF_dens, REF_pott, REF_qv, REF_pres,        & ! (in)
       0.0_RP, 0.0_RP, nd_order, nd_sfc_fact, nd_use_rs, & ! (in)
       DAMP_var(:,:,:,1), DAMP_var(:,:,:,2), DAMP_var(:,:,:,3), DAMP_var(:,:,:,4), DAMP_var(:,:,:,5), DAMP_var(:,:,:,6:6+QA-1), & ! (in)
       DAMP_alpha(:,:,:,1), DAMP_alpha(:,:,:,2), DAMP_alpha(:,:,:,3), DAMP_alpha(:,:,:,4), DAMP_alpha(:,:,:,5), & ! (in)
       DAMP_alpha(:,:,:,6:6+QA-1),                  & ! (in)
       divdmp_coef,                                 & ! (in)
       flag_fct_rho, flag_fct_momentum, flag_fct_t, & ! (in)
       flag_fct_along_stream,                       & ! (in)
       .false.,                                     & ! (in)
       1.0_DP, 1.0_DP, 1                            ) ! (in)

  message = "iq = ??"
  do iq = 1, QA
     write(message(6:7), "(i2)") iq
     write(message(9:11),"(a3)") "MAX"
     MINMAX(:,:,:) = Q_MAX * ( 1_RP + epsilon )
     call AssertLessThan(message, MINMAX(KS:KE,IS:IE,JS:JE), QTRC(KS:KE,IS:IE,JS:JE,iq))
     write(message(9:11),"(a3)") "MIN"
     MINMAX(:,:,:) = Q_MIN - epsilon
     call AssertGreaterThan(message, MINMAX(KS:KE,IS:IE,JS:JE), QTRC(KS:KE,IS:IE,JS:JE,iq))
  end do

end subroutine test_fctminmax


! private
subroutine copy
  DENS_o(:,:,:) = DENS(:,:,:)
  MOMZ_o(:,:,:) = MOMZ(:,:,:)
  MOMX_o(:,:,:) = MOMX(:,:,:)
  MOMY_o(:,:,:) = MOMY(:,:,:)
  RHOT_o(:,:,:) = RHOT(:,:,:)
  QTRC_o(:,:,:,:) = QTRC(:,:,:,:)
end subroutine copy


end module test_atmos_dyn
