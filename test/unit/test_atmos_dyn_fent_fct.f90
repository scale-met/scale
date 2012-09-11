module test_atmos_dyn_fent_fct

  !-----------------------------------------------------------------------------
  use mod_atmos_dyn, only: &
     ATMOS_DYN_main
  use dc_test, only: &
     AssertEqual, &
     AssertLessThan
    use mod_grid, only : &
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
  public :: test_atmos_dyn_fent_fct_run

  !-----------------------------------------------------------------------------
  include 'inc_index.h'
  include 'inc_tracer.h'
  include 'inc_precision.h'
  !-----------------------------------------------------------------------------
  real(RP) :: DENS(KA,IA,JA)
  real(RP) :: MOMZ(KA,IA,JA)
  real(RP) :: MOMX(KA,IA,JA)
  real(RP) :: MOMY(KA,IA,JA)
  real(RP) :: RHOT(KA,IA,JA)
  real(RP) :: QTRC(KA,IA,JA,QA)

  real(RP) :: DENS_o(KA,IA,JA)
  real(RP) :: MOMZ_o(KA,IA,JA)
  real(RP) :: MOMX_o(KA,IA,JA)
  real(RP) :: MOMY_o(KA,IA,JA)
  real(RP) :: RHOT_o(KA,IA,JA)
  real(RP) :: QTRC_o(KA,IA,JA,QA)

  real(RP) :: QDRY(KA,IA,JA)
  real(RP) :: DDIV(KA,IA,JA)
  real(RP) :: SINK(KA,IA,JA,5+QA)


  real(RP) :: REF_dens(KA)
  real(RP) :: REF_pott(KA)

  real(RP) :: DAMP_var(KA,IA,JA,5)
  real(RP) :: DAMP_alpha(KA,IA,JA,5)

  real(RP) :: CNDZ(3,KA)
  real(RP) :: CNMZ(3,KA)
  real(RP) :: CNDX(3,IA)
  real(RP) :: CNMX(3,IA)
  real(RP) :: CNDY(3,JA)
  real(RP) :: CNMY(3,JA)

  real(RP) :: CORIOLI(1,IA,JA)

  real(RP) :: DIFF4, DIFF2
  real(RP) :: divdmp_coef, LSsink_D

  real(RP), save :: ZERO(KA,IA,JA)

  integer, save :: KME ! end of main region

  integer :: k, i, j, iq
  character(len=7) :: message
  !-----------------------------------------------------------------------------
contains

  subroutine test_atmos_dyn_fent_fct_run
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_atmos_dyn, only: &
     ATMOS_DYN_init
  use mod_grid, only: &
     GRID_CZ_mask

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
  ZERO(:,:,:) = 0.0_RP

  do j = 1, JA
     lat(1,:,j) = real(j, RP)
  end do
  call ATMOS_DYN_init( DIFF4, DIFF2, CORIOLI, SINK,         & ! (out)
                       CNDZ, CNMZ, CNDX, CNMX, CNDY, CNMY,  & ! (out)
                       CDZ, CDX, CDY,                       & ! (in)
                       lat,                                 & ! (in)
                       0.0_RP, 0.0_RP, .false.              ) ! (in)


  do k = KS+1, KE
     if ( .not. GRID_CZ_mask(k) ) then
        KME = k - 1
        exit
     end if
  end do

  divdmp_coef = 0.0_RP
  LSsink_D    = 0.0_RP

  !########## test ##########


  call test_undef

  call test_conserve

end subroutine test_atmos_dyn_fent_fct_run
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

  do k = 1, KA
     REF_dens(k) = KA - k + 1
     REF_pott(k) = k
  end do

  DAMP_var  (:,:,:,:) = -9.999E30_RP
  DAMP_alpha(:,:,:,:) = 0.0_RP

  do i = 1, 2
     call ATMOS_DYN_main( &
          DENS, MOMZ, MOMX, MOMY, RHOT, QTRC,   & ! (inout)
          SINK,                                 & ! (inout)
          QDRY, DDIV,                           & ! (out)
          CNDZ, CNMZ, CNDX, CNMX, CNDY, CNMY,   & ! (in)
          CZ, FZ, CDZ, CDX, CDY, FDZ, FDX, FDY, & ! (in)
          RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,   & ! (in)
          REF_dens, REF_pott, DIFF4, DIFF2,     & ! (in)
          DAMP_var, DAMP_alpha,                 & ! (in)
          divdmp_coef, LSsink_D,                & ! (in)
          1.0_RP, 1                             ) ! (in)
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

subroutine test_conserve
  use mod_atmos_vars, only: &
       ATMOS_vars_fillhalo

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

  REF_dens(:) = 1.0_RP
  REF_pott(:) = 1.0_RP

  DAMP_var  (:,:,:,:) = -9.999E30_RP
  DAMP_alpha(:,:,:,:) = 0.0_RP

  call ATMOS_vars_fillhalo

  call copy

  call ATMOS_DYN_main( &
         DENS, MOMZ, MOMX, MOMY, RHOT, QTRC,   & ! (inout)
         SINK,                                 & ! (inout)
         QDRY, DDIV,                           & ! (out)
         CNDZ, CNMZ, CNDX, CNMX, CNDY, CNMY,   & ! (in)
         CZ, FZ, CDZ, CDX, CDY, FDZ, FDX, FDY, & ! (in)
         RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,   & ! (in)
         REF_dens, REF_pott, DIFF4, DIFF2,     & ! (in)
         DAMP_var, DAMP_alpha,                 & ! (in)
         divdmp_coef, LSsink_D,                & ! (in)
         1.0_RP, 1                             ) ! (in)

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
  call AssertEqual("DENS", total_o, total, RP*2-2, 10)

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
  call AssertEqual("MOMX", total_o, total, RP*2-2, 10)

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
  call AssertEqual("MOMY", total_o, total, RP*2-2, 10)

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
  call AssertEqual("RHOT", total_o, total, RP*2-2, 10)

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
     call AssertEqual(message, total_o, total, RP*2-2, 10)
  end do

end subroutine test_conserve



! private
subroutine copy
  DENS_o(:,:,:) = DENS(:,:,:)
  MOMZ_o(:,:,:) = MOMZ(:,:,:)
  MOMX_o(:,:,:) = MOMX(:,:,:)
  MOMY_o(:,:,:) = MOMY(:,:,:)
  RHOT_o(:,:,:) = RHOT(:,:,:)
  QTRC_o(:,:,:,:) = QTRC(:,:,:,:)
end subroutine copy


end module test_atmos_dyn_fent_fct
