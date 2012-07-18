module test_atmos_phy_tb_smg

  !-----------------------------------------------------------------------------
  use mod_atmos_phy_tb, only: &
     ATMOS_PHY_TB_main
  use dc_test, only: &
     AssertEqual, &
     AssertLessThan
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  public :: test_atmos_phy_tb_smg_run

  !-----------------------------------------------------------------------------
  include 'inc_index.h'
  include 'inc_tracer.h'
  include 'inc_precision.h'
  !-----------------------------------------------------------------------------
  real(RP) :: MOMZ_t(KA,IA,JA)
  real(RP) :: MOMX_t(KA,IA,JA)
  real(RP) :: MOMY_t(KA,IA,JA)
  real(RP) :: RHOT_t(KA,IA,JA)
  real(RP) :: QTRC_t(KA,IA,JA,QA)

  real(RP) :: tke (KA,IA,JA) ! TKE
  real(RP) :: nu_C(KA,IA,JA) ! eddy viscosity (center)
  real(RP) :: Pr  (KA,IA,JA) ! Prantle number
  real(RP) :: Ri  (KA,IA,JA) ! Richardson number

  real(RP) :: MOMZ(KA,IA,JA)
  real(RP) :: MOMX(KA,IA,JA)
  real(RP) :: MOMY(KA,IA,JA)
  real(RP) :: RHOT(KA,IA,JA)
  real(RP) :: DENS(KA,IA,JA)
  real(RP) :: QTRC(KA,IA,JA,QA)

  real(RP) :: SFLX_MOMZ(IA,JA)
  real(RP) :: SFLX_MOMX(IA,JA)
  real(RP) :: SFLX_MOMY(IA,JA)
  real(RP) :: SFLX_POTT(IA,JA)
  real(RP) :: SFLX_QV  (IA,JA)

  real(RP), save :: ZERO(KA,IA,JA)

  integer, save :: KME ! end of main region

  integer :: k, i, j, iq
  character(len=7) :: message
  !-----------------------------------------------------------------------------
contains

  subroutine test_atmos_phy_tb_smg_run
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_atmos_phy_tb, only: &
     ATMOS_PHY_TB_setup
  use mod_grid, only: &
     GRID_CZ_mask

  !-----------------------------------------------------------------------------
  implicit none
  !-----------------------------------------------------------------------------
  !
  !++ parameters & variables
  !
  !=============================================================================


  !########## Initial setup ##########
  call ATMOS_PHY_TB_setup

  ZERO(:,:,:) = 0.0_RP

  do k = KS+1, KE
     if ( .not. GRID_CZ_mask(k) ) then
        KME = k - 1
        exit
     end if
  end do

  !########## test ##########

  call test_zero

  call test_constant

  call test_linear_z
  call test_linear_x
  call test_linear_y

  call test_energy

  call test_big

  call test_double

  call test_sfc_flux

  call test_noise

end subroutine test_atmos_phy_tb_smg_run
!=============================================================================


subroutine test_zero

  write(*,*) "Test zero"

  SFLX_MOMZ(:,:) = 0.0_RP
  SFLX_MOMX(:,:) = 0.0_RP
  SFLX_MOMY(:,:) = 0.0_RP
  SFLX_POTT(:,:) = 0.0_RP
  SFLX_QV  (:,:) = 0.0_RP

  MOMZ(:,:,:) = 0.0_RP
  MOMX(:,:,:) = 0.0_RP
  MOMY(:,:,:) = 0.0_RP
  RHOT(:,:,:) = 1.0_RP
  DENS(:,:,:) = 1.0_RP
  QTRC(:,:,:,:) = 0.0_RP

  call ATMOS_PHY_TB_main( &
       MOMZ_t, MOMX_t, MOMY_t, RHOT_t, QTRC_t, & ! (out)
       tke, nu_C, Ri, Pr,                      & ! (out)
       MOMZ, MOMX, MOMY, RHOT, DENS, QTRC,     & ! (in)
       SFLX_MOMZ, SFLX_MOMX, SFLX_MOMY, SFLX_POTT, SFLX_QV & ! (in)
       )

  call AssertEqual("MOMZ_t", ZERO(KS:KE,IS:IE,JS:JE), MOMZ_t(KS:KE,IS:IE,JS:JE))
  call AssertEqual("MOMX_t", ZERO(KS:KE,IS:IE,JS:JE), MOMX_t(KS:KE,IS:IE,JS:JE))
  call AssertEqual("MOMY_t", ZERO(KS:KE,IS:IE,JS:JE), MOMY_t(KS:KE,IS:IE,JS:JE))
  call AssertEqual("RHOT_t", ZERO(KS:KE,IS:IE,JS:JE), RHOT_t(KS:KE,IS:IE,JS:JE))
  message = "iq = ??"
  do iq = 1, QA
     write(message(6:7), "(i2)") iq
     call AssertEqual(message, ZERO(KS:KE,IS:IE,JS:JE), RHOT_t(KS:KE,IS:IE,JS:JE))
  end do

end subroutine test_zero
!=============================================================================
subroutine test_constant

  write(*,*) "Test constant"

  SFLX_MOMZ(:,:) = 0.0_RP
  SFLX_MOMX(:,:) = 0.0_RP
  SFLX_MOMY(:,:) = 0.0_RP
  SFLX_POTT(:,:) = 0.0_RP
  SFLX_QV  (:,:) = 0.0_RP

  MOMZ(:,:,:) = 1.0_RP
  MOMX(:,:,:) = 1.0_RP
  MOMY(:,:,:) = 1.0_RP
  RHOT(:,:,:) = 1.0_RP
  DENS(:,:,:) = 1.0_RP
  QTRC(:,:,:,:) = 1.0_RP

  call fill_halo(MOMZ, MOMX, MOMY, RHOT, DENS, QTRC)

  call ATMOS_PHY_TB_main( &
       MOMZ_t, MOMX_t, MOMY_t, RHOT_t, QTRC_t, & ! (out)
       tke, nu_C, Ri, Pr,                      & ! (out)
       MOMZ, MOMX, MOMY, RHOT, DENS, QTRC,     & ! (in)
       SFLX_MOMZ, SFLX_MOMX, SFLX_MOMY, SFLX_POTT, SFLX_QV & ! (in)
       )

  call AssertEqual("MOMZ_t", ZERO(KS+1:KE,IS:IE,JS:JE), MOMZ_t(KS+1:KE,IS:IE,JS:JE))
  call AssertEqual("MOMX_t", ZERO(KS:KE,IS:IE,JS:JE), MOMX_t(KS:KE,IS:IE,JS:JE))
  call AssertEqual("MOMY_t", ZERO(KS:KE,IS:IE,JS:JE), MOMY_t(KS:KE,IS:IE,JS:JE))
  call AssertEqual("RHOT_t", ZERO(KS:KE,IS:IE,JS:JE), RHOT_t(KS:KE,IS:IE,JS:JE))
  message = "iq = ??"
  do iq = 1, QA
     write(message(6:7), "(i2)") iq
     call AssertEqual(message, ZERO(KS:KE,IS:IE,JS:JE), RHOT_t(KS:KE,IS:IE,JS:JE))
  end do

end subroutine test_constant
!=============================================================================
subroutine test_linear_z
  use mod_grid, only: &
       GRID_CZ
  implicit none

  write(*,*) "Test linear z"

  SFLX_MOMZ(:,:) = 0.0_RP
  SFLX_MOMX(:,:) = 0.0_RP
  SFLX_MOMY(:,:) = 0.0_RP
  SFLX_POTT(:,:) = 0.0_RP
  SFLX_QV  (:,:) = 0.0_RP

  do k = 1, KA
     MOMZ(k,:,:) = 2.0_RP * GRID_CZ(k)
  end do
  MOMX(:,:,:) = 0.0_RP
  MOMY(:,:,:) = 0.0_RP
  RHOT(:,:,:) = 1.0_RP
  DENS(:,:,:) = 1.0_RP
  QTRC(:,:,:,:) = 0.0_RP

  call fill_halo(MOMZ, MOMX, MOMY, RHOT, DENS, QTRC)

  call ATMOS_PHY_TB_main( &
       MOMZ_t, MOMX_t, MOMY_t, RHOT_t, QTRC_t, & ! (out)
       tke, nu_C, Ri, Pr,                      & ! (out)
       MOMZ, MOMX, MOMY, RHOT, DENS, QTRC,     & ! (in)
       SFLX_MOMZ, SFLX_MOMX, SFLX_MOMY, SFLX_POTT, SFLX_QV & ! (in)
       )

  call AssertEqual("MOMZ_t", ZERO(KS+1:KME-1,IS:IE,JS:JE), MOMZ_t(KS+1:KME-1,IS:IE,JS:JE))
  call AssertEqual("MOMX_t", ZERO(KS:KE,IS:IE,JS:JE), MOMX_t(KS:KE,IS:IE,JS:JE))
  call AssertEqual("MOMY_t", ZERO(KS:KE,IS:IE,JS:JE), MOMY_t(KS:KE,IS:IE,JS:JE))
  call AssertEqual("RHOT_t", ZERO(KS:KE,IS:IE,JS:JE), RHOT_t(KS:KE,IS:IE,JS:JE))
  message = "iq = ??"
  do iq = 1, QA
     write(message(6:7), "(i2)") iq
     call AssertEqual(message, ZERO(KS:KE,IS:IE,JS:JE), RHOT_t(KS:KE,IS:IE,JS:JE))
  end do

end subroutine test_linear_z
!=============================================================================
subroutine test_linear_x
  use mod_grid, only: &
       GRID_CX

  write(*,*) "Test linear x"

  SFLX_MOMZ(:,:) = 0.0_RP
  SFLX_MOMX(:,:) = 0.0_RP
  SFLX_MOMY(:,:) = 0.0_RP
  SFLX_POTT(:,:) = 0.0_RP
  SFLX_QV  (:,:) = 0.0_RP

  do i = 1, IA
     MOMX(:,i,:) = 2.0_RP * GRID_CX(i)
  end do
  MOMZ(:,:,:) = 0.0_RP
  MOMY(:,:,:) = 0.0_RP
  RHOT(:,:,:) = 1.0_RP
  DENS(:,:,:) = 1.0_RP
  QTRC(:,:,:,:) = 0.0_RP

  call fill_halo(MOMZ, MOMX, MOMY, RHOT, DENS, QTRC)

  call ATMOS_PHY_TB_main( &
       MOMZ_t, MOMX_t, MOMY_t, RHOT_t, QTRC_t, & ! (out)
       tke, nu_C, Ri, Pr,                      & ! (out)
       MOMZ, MOMX, MOMY, RHOT, DENS, QTRC,     & ! (in)
       SFLX_MOMZ, SFLX_MOMX, SFLX_MOMY, SFLX_POTT, SFLX_QV & ! (in)
       )

  call AssertEqual("MOMZ_t", ZERO(KS+1:KME-1,IS:IE,JS:JE), MOMZ_t(KS+1:KME-1,IS:IE,JS:JE))
  call AssertEqual("MOMX_t", ZERO(KS:KE,IS:IE,JS:JE), MOMX_t(KS:KE,IS:IE,JS:JE))
  call AssertEqual("MOMY_t", ZERO(KS:KE,IS:IE,JS:JE), MOMY_t(KS:KE,IS:IE,JS:JE))
  call AssertEqual("RHOT_t", ZERO(KS:KE,IS:IE,JS:JE), RHOT_t(KS:KE,IS:IE,JS:JE))
  message = "iq = ??"
  do iq = 1, QA
     write(message(6:7), "(i2)") iq
     call AssertEqual(message, ZERO(KS:KE,IS:IE,JS:JE), RHOT_t(KS:KE,IS:IE,JS:JE))
  end do

end subroutine test_linear_x
!=============================================================================
subroutine test_linear_y
  use mod_grid, only: &
       GRID_CY

  write(*,*) "Test linear y"

  SFLX_MOMZ(:,:) = 0.0_RP
  SFLX_MOMX(:,:) = 0.0_RP
  SFLX_MOMY(:,:) = 0.0_RP
  SFLX_POTT(:,:) = 0.0_RP
  SFLX_QV  (:,:) = 0.0_RP

  do j = 1, JA
     MOMY(:,:,j) = 2.0_RP * GRID_CY(j)
  end do
  MOMZ(:,:,:) = 0.0_RP
  MOMX(:,:,:) = 0.0_RP
  RHOT(:,:,:) = 1.0_RP
  DENS(:,:,:) = 1.0_RP
  QTRC(:,:,:,:) = 0.0_RP

  call fill_halo(MOMZ, MOMX, MOMY, RHOT, DENS, QTRC)

  call ATMOS_PHY_TB_main( &
       MOMZ_t, MOMX_t, MOMY_t, RHOT_t, QTRC_t, & ! (out)
       tke, nu_C, Ri, Pr,                      & ! (out)
       MOMZ, MOMX, MOMY, RHOT, DENS, QTRC,     & ! (in)
       SFLX_MOMZ, SFLX_MOMX, SFLX_MOMY, SFLX_POTT, SFLX_QV & ! (in)
       )

  call AssertEqual("MOMZ_t", ZERO(KS+1:KME-1,IS:IE,JS:JE), MOMZ_t(KS+1:KME-1,IS:IE,JS:JE))
  call AssertEqual("MOMX_t", ZERO(KS:KE,IS:IE,JS:JE), MOMX_t(KS:KE,IS:IE,JS:JE))
  call AssertEqual("MOMY_t", ZERO(KS:KE,IS:IE,JS:JE), MOMY_t(KS:KE,IS:IE,JS:JE))
  call AssertEqual("RHOT_t", ZERO(KS:KE,IS:IE,JS:JE), RHOT_t(KS:KE,IS:IE,JS:JE))
  message = "iq = ??"
  do iq = 1, QA
     write(message(6:7), "(i2)") iq
     call AssertEqual(message, ZERO(KS:KE,IS:IE,JS:JE), RHOT_t(KS:KE,IS:IE,JS:JE))
  end do

end subroutine test_linear_y
!=============================================================================
subroutine test_energy
  use mod_grid, only: &
       GRID_FZ, &
       GRID_FX, &
       GRID_FY

  real(RP) :: eng1, eng2
  real(RP) :: PI2
  real(RP) :: dt
  PI2 = atan( 1.0_RP )*8.0_RP

  write(*,*) "Test energy"
  ! Smagorinsky-Lilly SGS model has no back scatter,
  ! whcih means that energy of resolved grid scale must decrease.

  SFLX_MOMZ(:,:) = 0.0_RP
  SFLX_MOMX(:,:) = 0.0_RP
  SFLX_MOMY(:,:) = 0.0_RP
  SFLX_POTT(:,:) = 0.0_RP
  SFLX_QV  (:,:) = 0.0_RP

  do j = 1, JA
  do i = 1, IA
  do k = 1, KA
     MOMZ(k,i,j) = 1.0_RP * sin( PI2 * ( k*1.0_RP/(KE-KS+1) + i*2.0_RP/(IE-IS+1) + j*3.0_RP/(JE-JS+1) ) )
     MOMX(k,i,j) = 2.0_RP * cos( PI2 * ( k*2.0_RP/(KE-KS+1) + i*3.0_RP/(IE-IS+1) + j*1.0_RP/(JE-JS+1) ) )
     MOMY(k,i,j) = 3.0_RP * sin( PI2 * ( k*3.0_RP/(KE-KS+1) + i*1.0_RP/(IE-IS+1) + j*2.0_RP/(JE-JS+1) ) )
  end do
  end do
  end do
  RHOT(:,:,:) = 1.0_RP
  DENS(:,:,:) = 1.0_RP
  QTRC(:,:,:,:) = 0.0_RP

  dt = 1.0E-3_RP *  min( min( GRID_FZ(KE)-GRID_FZ(KS-1), GRID_FX(IE)-GRID_FX(IS-1) ), GRID_FY(JE)-GRID_FY(JS-1) ) / ( PI2 * 3.0_RP )

  eng1 = 0.0_RP
  do j = JS, JE
  do i = IS, IE
  do k = KS, KE
     eng1 = eng1 + &
          0.5_RP * ( (MOMZ(k,i,j)/DENS(k,i,j))**2 + (MOMX(k,i,j)/DENS(k,i,j))**2 + (MOMY(k,i,j)/DENS(k,i,j))**2 )
  end do
  end do
  end do

  call fill_halo(MOMZ, MOMX, MOMY, RHOT, DENS, QTRC)

  call ATMOS_PHY_TB_main( &
       MOMZ_t, MOMX_t, MOMY_t, RHOT_t, QTRC_t, & ! (out)
       tke, nu_C, Ri, Pr,                      & ! (out)
       MOMZ, MOMX, MOMY, RHOT, DENS, QTRC,     & ! (in)
       SFLX_MOMZ, SFLX_MOMX, SFLX_MOMY, SFLX_POTT, SFLX_QV & ! (in)
       )

  do j = JS, JE
  do i = IS, IE
  do k = KS, KE
     MOMZ(k,i,j) = MOMZ(k,i,j) + dt * MOMZ_t(k,i,j)
     MOMX(k,i,j) = MOMX(k,i,j) + dt * MOMX_t(k,i,j)
     MOMY(k,i,j) = MOMY(k,i,j) + dt * MOMY_t(k,i,j)
  end do
  end do
  end do

  eng2 = 0.0_RP
  do j = JS, JE
  do i = IS, IE
  do k = KS, KE
     eng2 = eng2 + &
          0.5_RP * ( (MOMZ(k,i,j)/DENS(k,i,j))**2 + (MOMX(k,i,j)/DENS(k,i,j))**2 + (MOMY(k,i,j)/DENS(k,i,j))**2 )
  end do
  end do
  end do

  call AssertLessThan("energy", eng1, eng2)

end subroutine test_energy
!=============================================================================
subroutine test_big
  use mod_grid, only: &
       GRID_FZ, &
       GRID_FX, &
       GRID_FY

  real(RP) :: BIG(KA,IA,JA)

  real(RP) :: eng1, eng2
  real(RP) :: PI2
  real(RP) :: dt
  PI2 = atan( 1.0_RP )*8.0_RP

  BIG(:,:,:) = 9.99E8_RP

  write(*,*) "Test big"
  ! check not to include BUG (UNDEF) value

  SFLX_MOMZ(:,:) = 0.0_RP
  SFLX_MOMX(:,:) = 0.0_RP
  SFLX_MOMY(:,:) = 0.0_RP
  SFLX_POTT(:,:) = 0.0_RP
  SFLX_QV  (:,:) = 0.0_RP

  do j = 1, JA
  do i = 1, IA
  do k = 1, KA
     MOMZ(k,i,j) = 1.0_RP * sin( k*1.0_RP + i*2.0_RP + j*3.0_RP )
     MOMX(k,i,j) = 2.0_RP * cos( k*2.0_RP + i*3.0_RP + j*1.0_RP )
     MOMY(k,i,j) = 3.0_RP * sin( k*3.0_RP + i*1.0_RP + j*2.0_RP )
     RHOT(k,i,j) = 4.0_RP * cos( k*1.0_RP + i*1.0_RP + j*3.0_RP ) + 300.0_RP
     DENS(k,i,j) = 5.0_RP * sin( k*2.0_RP + i*2.0_RP + j*2.0_RP ) + 6.0_rp
     do iq = 1, QA
        QTRC(k,i,j,iq) = 6.0_RP * sin( k*3.0_RP + i*3.0_RP + j*1.0_RP + iq*2.0_RP ) + 6.0_RP
     end do
  end do
  end do
  end do

  call fill_halo(MOMZ, MOMX, MOMY, RHOT, DENS, QTRC)

  call ATMOS_PHY_TB_main( &
       MOMZ_t, MOMX_t, MOMY_t, RHOT_t, QTRC_t, & ! (out)
       tke, nu_C, Ri, Pr,                      & ! (out)
       MOMZ, MOMX, MOMY, RHOT, DENS, QTRC,     & ! (in)
       SFLX_MOMZ, SFLX_MOMX, SFLX_MOMY, SFLX_POTT, SFLX_QV & ! (in)
       )

  call AssertLessThan("MOMZ_t", BIG(KS:KE,IS:IE,JS:JE), abs(MOMZ_t(KS:KE,IS:IE,JS:JE)))
  call AssertLessThan("MOMX_t", BIG(KS:KE,IS:IE,JS:JE), abs(MOMX_t(KS:KE,IS:IE,JS:JE)))
  call AssertLessThan("MOMY_t", BIG(KS:KE,IS:IE,JS:JE), abs(MOMY_t(KS:KE,IS:IE,JS:JE)))
  call AssertLessThan("RHOT_t", BIG(KS:KE,IS:IE,JS:JE), abs(RHOT_t(KS:KE,IS:IE,JS:JE)))
  message = "iq = ??"
  do iq = 1, QA
     write(message(6:7), "(i2)") iq
     call AssertLessThan(message, BIG(KS:KE,IS:IE,JS:JE), abs(QTRC_t(KS:KE,IS:IE,JS:JE,iq)))
  end do

end subroutine test_big
!=============================================================================
subroutine test_noise

  write(*,*) "Test Noise"
  ! Smagorinsky-Lilly SGS model has no back scatter,
  ! whcih means that two-grid noise must be reduced

  SFLX_MOMZ(:,:) = 0.0_RP
  SFLX_MOMX(:,:) = 0.0_RP
  SFLX_MOMY(:,:) = 0.0_RP
  SFLX_POTT(:,:) = 0.0_RP
  SFLX_QV  (:,:) = 0.0_RP

  do j = 1, JA
  do i = 1, IA
  do k = KS, KE
     MOMZ(k,i,j) = real( mod(k+i+j,2) * 2 - 1, RP ) ! -1 or 1
     MOMX(k,i,j) = MOMZ(k,i,j)
     MOMY(k,i,j) = MOMZ(k,i,j)
     RHOT(k,i,j) = MOMZ(k,i,j) + 300.0_RP
     QTRC(k,i,j,:) = MOMZ(k,i,j)
  end do
  end do
  end do
  DENS(:,:,:) = 1.0_RP

  call fill_halo(MOMZ, MOMX, MOMY, RHOT, DENS, QTRC)

  call ATMOS_PHY_TB_main( &
       MOMZ_t, MOMX_t, MOMY_t, RHOT_t, QTRC_t, & ! (out)
       tke, nu_C, Ri, Pr,                      & ! (out)
       MOMZ, MOMX, MOMY, RHOT, DENS, QTRC,     & ! (in)
       SFLX_MOMZ, SFLX_MOMX, SFLX_MOMY, SFLX_POTT, SFLX_QV & ! (in)
       )

  ! tendency must be opposit sign
  call AssertLessThan("MOMZ", ZERO(KS:KE,IS:IE,JS:JE)+1.E-10_RP, MOMZ_t(KS:KE,IS:IE,JS:JE)*MOMZ(KS:KE,IS:IE,JS:JE) )
  call AssertLessThan("MOMX", ZERO(KS:KE,IS:IE,JS:JE)+1.E-10_RP, MOMX_t(KS:KE,IS:IE,JS:JE)*MOMX(KS:KE,IS:IE,JS:JE) )
  call AssertLessThan("MOMY", ZERO(KS:KE,IS:IE,JS:JE)+1.E-10_RP, MOMY_t(KS:KE,IS:IE,JS:JE)*MOMY(KS:KE,IS:IE,JS:JE) )
  call AssertLessThan("RHOT", ZERO(KS:KE,IS:IE,JS:JE)+1.E-10_RP, RHOT_t(KS:KE,IS:IE,JS:JE)*(RHOT(KS:KE,IS:IE,JS:JE)-300.0_RP) )
  message = "iq = ??"
  do iq = 1, QA
     write(message(6:7), "(i2)") iq
     call AssertLessThan(message, ZERO(KS:KE,IS:IE,JS:JE)+1.E-10_RP, QTRC_t(KS:KE,IS:IE,JS:JE,iq)*QTRC(KS:KE,IS:IE,JS:JE,iq) )
  end do


end subroutine test_noise
!=============================================================================
subroutine test_double
  use mod_grid, only: &
       GRID_FZ, &
       GRID_FX, &
       GRID_FY

  real(RP) :: MOMZ_t2(KA,IA,JA)
  real(RP) :: MOMX_t2(KA,IA,JA)
  real(RP) :: MOMY_t2(KA,IA,JA)
  real(RP) :: RHOT_t2(KA,IA,JA)
  real(RP) :: QTRC_t2(KA,IA,JA,QA)

  real(RP) :: FOUR(KA,IA,JA)
  real(RP) :: PI2

  FOUR(:,:,:) = 4.0_RP
  PI2 = atan( 1.0_RP )*8.0_RP

  write(*,*) "Test double"

  SFLX_MOMZ(:,:) = 0.0_RP
  SFLX_MOMX(:,:) = 0.0_RP
  SFLX_MOMY(:,:) = 0.0_RP
  SFLX_POTT(:,:) = 0.0_RP
  SFLX_QV  (:,:) = 0.0_RP

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

  call fill_halo(MOMZ, MOMX, MOMY, RHOT, DENS, QTRC)

  call ATMOS_PHY_TB_main( &
       MOMZ_t, MOMX_t, MOMY_t, RHOT_t, QTRC_t, & ! (out)
       tke, nu_C, Ri, Pr,                      & ! (out)
       MOMZ, MOMX, MOMY, RHOT, DENS, QTRC,     & ! (in)
       SFLX_MOMZ, SFLX_MOMX, SFLX_MOMY, SFLX_POTT, SFLX_QV & ! (in)
       )

  MOMZ(:,:,:) = MOMZ(:,:,:) * 2.0_RP
  MOMX(:,:,:) = MOMX(:,:,:) * 2.0_RP
  MOMY(:,:,:) = MOMY(:,:,:) * 2.0_RP
  QTRC(:,:,:,:) = QTRC(:,:,:,:) * 2.0_RP

  call ATMOS_PHY_TB_main( &
       MOMZ_t2, MOMX_t2, MOMY_t2, RHOT_t2, QTRC_t2, & ! (out)
       tke, nu_C, Ri, Pr,                      & ! (out)
       MOMZ, MOMX, MOMY, RHOT, DENS, QTRC,     & ! (in)
       SFLX_MOMZ, SFLX_MOMX, SFLX_MOMY, SFLX_POTT, SFLX_QV & ! (in)
       )

  call AssertEqual("MOMZ_t", FOUR(KS+1:KME-1,IS:IE,JS:JE), MOMZ_t2(KS+1:KME-1,IS:IE,JS:JE)/MOMZ_t(KS+1:KME-1,IS:IE,JS:JE))
  call AssertEqual("MOMX_t", FOUR(KS:KE,IS:IE,JS:JE), MOMX_t2(KS:KE,IS:IE,JS:JE)/MOMX_t(KS:KE,IS:IE,JS:JE))
  call AssertEqual("MOMY_t", FOUR(KS:KE,IS:IE,JS:JE), MOMY_t2(KS:KE,IS:IE,JS:JE)/MOMY_t(KS:KE,IS:IE,JS:JE))
  message = "iq = ??"
  do iq = 1, QA
     write(message(6:7), "(i2)") iq
     call AssertEqual(message, FOUR(KS:KE,IS:IE,JS:JE), QTRC_t2(KS:KE,IS:IE,JS:JE,iq)/QTRC_t(KS:KE,IS:IE,JS:JE,iq))
  end do

end subroutine test_double
!=============================================================================
subroutine test_sfc_flux
  use mod_grid, only: &
       GRID_RFDZ

  write(*,*) "Test surface flux"

  do j = 1, JA
  do i = 1, IA
     SFLX_MOMZ(i,j) = real(i + j, RP) *     1.0_RP
     SFLX_MOMX(:,:) = real(i + j, RP) *    10.0_RP
     SFLX_MOMY(:,:) = real(i + j, RP) *   100.0_RP
     SFLX_POTT(:,:) = real(i + j, RP) *  1000.0_RP
     SFLX_QV  (:,:) = real(i + j, RP) * 10000.0_RP
  end do
  end do

  MOMZ(:,:,:) = 0.0_RP
  MOMX(:,:,:) = 0.0_RP
  MOMY(:,:,:) = 0.0_RP
  RHOT(:,:,:) = 1.0_RP
  DENS(:,:,:) = 1.0_RP
  QTRC(:,:,:,:) = 0.0_RP

  call ATMOS_PHY_TB_main( &
       MOMZ_t, MOMX_t, MOMY_t, RHOT_t, QTRC_t, & ! (out)
       tke, nu_C, Ri, Pr,                      & ! (out)
       MOMZ, MOMX, MOMY, RHOT, DENS, QTRC,     & ! (in)
       SFLX_MOMZ, SFLX_MOMX, SFLX_MOMY, SFLX_POTT, SFLX_QV & ! (in)
       )

  call AssertEqual("MOMZ_t", SFLX_MOMZ(IS:IE,JS:JE)*GRID_RFDZ(KS), MOMZ_t(KS,IS:IE,JS:JE) )
  call AssertEqual("MOMX_t", SFLX_MOMX(IS:IE,JS:JE)*GRID_RFDZ(KS), MOMX_t(KS,IS:IE,JS:JE) )
  call AssertEqual("MOMY_t", SFLX_MOMY(IS:IE,JS:JE)*GRID_RFDZ(KS), MOMY_t(KS,IS:IE,JS:JE) )
  call AssertEqual("POTT_t", SFLX_POTT(IS:IE,JS:JE)*GRID_RFDZ(KS), RHOT_t(KS,IS:IE,JS:JE)/DENS(KS,IS:IE,JS:JE) )
  call AssertEqual("QV",     SFLX_QV(IS:IE,JS:JE)*GRID_RFDZ(KS),   QTRC_t(KS,IS:IE,JS:JE,I_QV)/DENS(KS,IS:IE,JS:JE) )

end subroutine test_sfc_flux



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
