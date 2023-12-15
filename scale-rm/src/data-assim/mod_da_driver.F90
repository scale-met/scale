!-------------------------------------------------------------------------------
!> module Data Assimilation driver
!!
!! @par Description
!!          Data Assimilation driver
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_da_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: DA_driver_setup
  public :: DA_driver_finalize
  public :: DA_driver_update

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: make_ens_2D
  private :: make_ens_3D
  private :: make_output

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer :: datatype

  real(RP), allocatable :: QHYD(:,:,:,:)

  real(RP), allocatable :: RH   (:,:,:)
  real(RP), allocatable :: Qdry (:,:,:)
  real(RP), allocatable :: Rtot (:,:,:)
  real(RP), allocatable :: CVtot(:,:,:)
  real(RP), allocatable :: CPtot(:,:,:)

  real(RP), allocatable :: PREC(:,:)

  real(RP), allocatable :: ENS_mean_DENS(:,:,:)
  real(RP), allocatable :: ENS_mean_MOMX(:,:,:)
  real(RP), allocatable :: ENS_mean_MOMY(:,:,:)
  real(RP), allocatable :: ENS_mean_MOMZ(:,:,:)
  real(RP), allocatable :: ENS_mean_RHOT(:,:,:)
  real(RP), allocatable :: ENS_mean_U   (:,:,:)
  real(RP), allocatable :: ENS_mean_V   (:,:,:)
  real(RP), allocatable :: ENS_mean_W   (:,:,:)
  real(RP), allocatable :: ENS_mean_TEMP(:,:,:)
  real(RP), allocatable :: ENS_mean_PRES(:,:,:)
  real(RP), allocatable :: ENS_mean_QV  (:,:,:)
  real(RP), allocatable :: ENS_mean_QC  (:,:,:)
  real(RP), allocatable :: ENS_mean_QR  (:,:,:)
  real(RP), allocatable :: ENS_mean_QI  (:,:,:)
  real(RP), allocatable :: ENS_mean_QS  (:,:,:)
  real(RP), allocatable :: ENS_mean_QG  (:,:,:)
  real(RP), allocatable :: ENS_mean_RH  (:,:,:)

  real(RP), allocatable :: ENS_sprd_DENS(:,:,:)
  real(RP), allocatable :: ENS_sprd_MOMX(:,:,:)
  real(RP), allocatable :: ENS_sprd_MOMY(:,:,:)
  real(RP), allocatable :: ENS_sprd_MOMZ(:,:,:)
  real(RP), allocatable :: ENS_sprd_RHOT(:,:,:)
  real(RP), allocatable :: ENS_sprd_U   (:,:,:)
  real(RP), allocatable :: ENS_sprd_V   (:,:,:)
  real(RP), allocatable :: ENS_sprd_W   (:,:,:)
  real(RP), allocatable :: ENS_sprd_TEMP(:,:,:)
  real(RP), allocatable :: ENS_sprd_PRES(:,:,:)
  real(RP), allocatable :: ENS_sprd_QV  (:,:,:)
  real(RP), allocatable :: ENS_sprd_QC  (:,:,:)
  real(RP), allocatable :: ENS_sprd_QR  (:,:,:)
  real(RP), allocatable :: ENS_sprd_QI  (:,:,:)
  real(RP), allocatable :: ENS_sprd_QS  (:,:,:)
  real(RP), allocatable :: ENS_sprd_QG  (:,:,:)
  real(RP), allocatable :: ENS_sprd_RH  (:,:,:)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine DA_driver_setup
    use mpi
    use scale_prc, only: &
       PRC_LOCAL_COMM_WORLD, &
       PRC_nprocs,           &
       PRC_myrank,           &
       PRC_abort
    use scale_prc_cartesC, only: &
       PRC_NUM_X, &
       PRC_NUM_Y
    use scale_atmos_grid_cartesC, only: &
       DX, &
       DY
    use scale_atmos_grid_cartesC_index, only: &
       KMAX,       &
       IMAX,       &
       JMAX,       &
       KA, KS, KE, &
       IA, IS, IE, &
       JA, JS, JE, &
       KHALO,      &
       IHALO,      &
       JHALO
    use scale_atmos_hydrometeor, only: &
       N_HYD
    use scale_topography, only: &
       TOPOGRAPHY_Zsfc
    use scale_comm_ensemble, only: &
       ENSEMBLE_world  => COMM_ENSEMBLE_world,  &
       ENSEMBLE_nprocs => COMM_ENSEMBLE_nprocs, &
       ENSEMBLE_myrank => COMM_ENSEMBLE_myrank, &
       COMM_ENSEMBLE_setup
    use scale_letkf, only: &
       LETKF_setup
    use mod_da_vars, only: &
       OBS_IN_NUM
    use mod_da_param_estimation, only: &
       DA_param_estimation_setup
    implicit none
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("DA_driver_setup",*) 'Setup'

    if     ( RP == SP ) then
       datatype = MPI_REAL
    else if( RP == DP ) then
       datatype = MPI_DOUBLE_PRECISION
    else
       LOG_ERROR("DA_driver_setup",*) 'The precision has not been implemented yet:', RP
       call PRC_abort
    endif

    call COMM_ENSEMBLE_setup

    LOG_INFO("DA_driver_setup",*) 'ENSEMBLE myrank/nprocs:', ENSEMBLE_myrank, ENSEMBLE_nprocs

    allocate( QHYD( KA, IA, JA, N_HYD ) )
    QHYD(:,:,:,:) = 0.0_RP

    allocate( RH   ( KA, IA, JA ) )
    allocate( Qdry ( KA, IA, JA ) )
    allocate( Rtot ( KA, IA, JA ) )
    allocate( CVtot( KA, IA, JA ) )
    allocate( CPtot( KA, IA, JA ) )
    RH   (:,:,:) = 0.0_RP
    Qdry (:,:,:) = 0.0_RP
    Rtot (:,:,:) = 0.0_RP
    CVtot(:,:,:) = 0.0_RP
    CPtot(:,:,:) = 0.0_RP

    allocate( PREC( IA, JA ) )
    PREC(:,:) = 0.0_RP

    allocate( ENS_mean_DENS(KA,IA,JA) )
    allocate( ENS_mean_MOMX(KA,IA,JA) )
    allocate( ENS_mean_MOMY(KA,IA,JA) )
    allocate( ENS_mean_MOMZ(KA,IA,JA) )
    allocate( ENS_mean_RHOT(KA,IA,JA) )
    allocate( ENS_mean_U   (KA,IA,JA) )
    allocate( ENS_mean_V   (KA,IA,JA) )
    allocate( ENS_mean_W   (KA,IA,JA) )
    allocate( ENS_mean_TEMP(KA,IA,JA) )
    allocate( ENS_mean_PRES(KA,IA,JA) )
    allocate( ENS_mean_QV  (KA,IA,JA) )
    allocate( ENS_mean_QC  (KA,IA,JA) )
    allocate( ENS_mean_QR  (KA,IA,JA) )
    allocate( ENS_mean_QI  (KA,IA,JA) )
    allocate( ENS_mean_QS  (KA,IA,JA) )
    allocate( ENS_mean_QG  (KA,IA,JA) )
    allocate( ENS_mean_RH  (KA,IA,JA) )
    ENS_mean_DENS(:,:,:) = 0.0_RP
    ENS_mean_MOMX(:,:,:) = 0.0_RP
    ENS_mean_MOMY(:,:,:) = 0.0_RP
    ENS_mean_MOMZ(:,:,:) = 0.0_RP
    ENS_mean_RHOT(:,:,:) = 0.0_RP
    ENS_mean_U   (:,:,:) = 0.0_RP
    ENS_mean_V   (:,:,:) = 0.0_RP
    ENS_mean_W   (:,:,:) = 0.0_RP
    ENS_mean_TEMP(:,:,:) = 0.0_RP
    ENS_mean_PRES(:,:,:) = 0.0_RP
    ENS_mean_QV  (:,:,:) = 0.0_RP
    ENS_mean_QC  (:,:,:) = 0.0_RP
    ENS_mean_QR  (:,:,:) = 0.0_RP
    ENS_mean_QI  (:,:,:) = 0.0_RP
    ENS_mean_QS  (:,:,:) = 0.0_RP
    ENS_mean_QG  (:,:,:) = 0.0_RP
    ENS_mean_RH  (:,:,:) = 0.0_RP

    allocate( ENS_sprd_DENS(KA,IA,JA) )
    allocate( ENS_sprd_MOMX(KA,IA,JA) )
    allocate( ENS_sprd_MOMY(KA,IA,JA) )
    allocate( ENS_sprd_MOMZ(KA,IA,JA) )
    allocate( ENS_sprd_RHOT(KA,IA,JA) )
    allocate( ENS_sprd_U   (KA,IA,JA) )
    allocate( ENS_sprd_V   (KA,IA,JA) )
    allocate( ENS_sprd_W   (KA,IA,JA) )
    allocate( ENS_sprd_TEMP(KA,IA,JA) )
    allocate( ENS_sprd_PRES(KA,IA,JA) )
    allocate( ENS_sprd_QV  (KA,IA,JA) )
    allocate( ENS_sprd_QC  (KA,IA,JA) )
    allocate( ENS_sprd_QR  (KA,IA,JA) )
    allocate( ENS_sprd_QI  (KA,IA,JA) )
    allocate( ENS_sprd_QS  (KA,IA,JA) )
    allocate( ENS_sprd_QG  (KA,IA,JA) )
    allocate( ENS_sprd_RH  (KA,IA,JA) )
    ENS_sprd_DENS(:,:,:) = 0.0_RP
    ENS_sprd_MOMX(:,:,:) = 0.0_RP
    ENS_sprd_MOMY(:,:,:) = 0.0_RP
    ENS_sprd_MOMZ(:,:,:) = 0.0_RP
    ENS_sprd_RHOT(:,:,:) = 0.0_RP
    ENS_sprd_U   (:,:,:) = 0.0_RP
    ENS_sprd_V   (:,:,:) = 0.0_RP
    ENS_sprd_W   (:,:,:) = 0.0_RP
    ENS_sprd_TEMP(:,:,:) = 0.0_RP
    ENS_sprd_PRES(:,:,:) = 0.0_RP
    ENS_sprd_QV  (:,:,:) = 0.0_RP
    ENS_sprd_QC  (:,:,:) = 0.0_RP
    ENS_sprd_QR  (:,:,:) = 0.0_RP
    ENS_sprd_QI  (:,:,:) = 0.0_RP
    ENS_sprd_QS  (:,:,:) = 0.0_RP
    ENS_sprd_QG  (:,:,:) = 0.0_RP
    ENS_sprd_RH  (:,:,:) = 0.0_RP

    call LETKF_setup( OBS_IN_NUM,           &
                      ENSEMBLE_world,       &
                      ENSEMBLE_nprocs,      &
                      ENSEMBLE_myrank,      &
                      PRC_LOCAL_COMM_WORLD, &
                      PRC_nprocs,           &
                      PRC_myrank,           &
                      PRC_NUM_X,            &
                      PRC_NUM_Y,            &
                      KA, KS, KE,           &
                      IA, IS, IE,           &
                      JA, JS, JE,           &
                      KMAX,                 &
                      IMAX,                 &
                      JMAX,                 &
                      KHALO,                &
                      IHALO,                &
                      JHALO,                &
                      DX,                   &
                      DY,                   &
                      TOPOGRAPHY_Zsfc(:,:)  )

    call DA_param_estimation_setup

    return
  end subroutine DA_driver_setup

  !-----------------------------------------------------------------------------
  !> finalize
  subroutine DA_driver_finalize
    use scale_comm_ensemble, only: &
      COMM_ENSEMBLE_finalize
    use scale_letkf, only: &
       LETKF_finalize
    use mod_da_param_estimation, only: &
       DA_param_estimation_finalize
    implicit none
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("DA_driver_finalize",*) 'Finalize'

    deallocate( QHYD )

    deallocate( RH    )
    deallocate( Qdry  )
    deallocate( Rtot  )
    deallocate( CVtot )
    deallocate( CPtot )

    deallocate( PREC )

    deallocate( ENS_mean_DENS )
    deallocate( ENS_mean_MOMX )
    deallocate( ENS_mean_MOMY )
    deallocate( ENS_mean_MOMZ )
    deallocate( ENS_mean_RHOT )
    deallocate( ENS_mean_U    )
    deallocate( ENS_mean_V    )
    deallocate( ENS_mean_W    )
    deallocate( ENS_mean_TEMP )
    deallocate( ENS_mean_PRES )
    deallocate( ENS_mean_QV   )
    deallocate( ENS_mean_QC   )
    deallocate( ENS_mean_QR   )
    deallocate( ENS_mean_QI   )
    deallocate( ENS_mean_QS   )
    deallocate( ENS_mean_QG   )
    deallocate( ENS_mean_RH   )

    deallocate( ENS_sprd_U    )
    deallocate( ENS_sprd_V    )
    deallocate( ENS_sprd_W    )
    deallocate( ENS_sprd_TEMP )
    deallocate( ENS_sprd_PRES )
    deallocate( ENS_sprd_QV   )
    deallocate( ENS_sprd_QC   )
    deallocate( ENS_sprd_QR   )
    deallocate( ENS_sprd_QI   )
    deallocate( ENS_sprd_QS   )
    deallocate( ENS_sprd_QG   )
    deallocate( ENS_sprd_RH   )

    call DA_param_estimation_finalize
    call LETKF_finalize
    call COMM_ENSEMBLE_finalize

    return
  end subroutine DA_driver_finalize

  !-----------------------------------------------------------------------------
  !> Data Assimilation
  subroutine DA_driver_update
    use scale_const, &
       PRE00 => CONST_PRE00
    use scale_atmos_grid_cartesC_index, only: &
       KA, KS, KE, &
       IA, IS, IE, &
       JA, JS, JE
    use scale_atmos_grid_cartesC_real, only: &
       HGT => ATMOS_GRID_CARTESC_REAL_CZ
    use scale_atmos_hydrometeor, only: &
       I_HC, &
       I_HR, &
       I_HI, &
       I_HS, &
       I_HG
    use scale_tracer, only: &
       QA,          &
       TRACER_CV,   &
       TRACER_CP,   &
       TRACER_R,    &
       TRACER_MASS
    use scale_atmos_thermodyn, only: &
       ATMOS_THERMODYN_specific_heat
    use scale_comm_ensemble, only: &
       ENSEMBLE_world  => COMM_ENSEMBLE_world,  &
       ENSEMBLE_nprocs => COMM_ENSEMBLE_nprocs, &
       ENSEMBLE_myrank => COMM_ENSEMBLE_myrank
    use scale_topography, only: &
       TOPO => TOPOGRAPHY_Zsfc
    use scale_letkf, only: &
       LETKF_obs_readfile,   &
       LETKF_obs_clear,      &
       LETKF_obs_operator,   &
       LETKF_obs_initialize, &
       LETKF_system
    use mod_atmos_vars, only: &
       ATMOS_vars_get_diagnostic, &
       U,    &
       V,    &
       W,    &
       TEMP, &
       PRES, &
       POTT, &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       RHOT, &
       QTRC, &
       QV,   &
       QC,   &
       QR,   &
       QI,   &
       QS,   &
       QG
    use mod_atmos_phy_sf_vars, only: &
       PS   => ATMOS_PHY_SF_SFC_PRES, &
       U10M => ATMOS_PHY_SF_U10,      &
       V10M => ATMOS_PHY_SF_V10,      &
       T2M  => ATMOS_PHY_SF_T2,       &
       Q2M  => ATMOS_PHY_SF_Q2
    use mod_atmos_phy_mp_driver, only: &
       ATMOS_PHY_MP_driver_qhyd2qtrc
    use mod_da_vars, only: &
       OBS_IN_NUM,             &
       OBS_IN_FORMAT,          &
       OBS_IN_BASENAME,        &
       OBS_IN_MASKFILE,        &
       POSITIVE_DEFINITE_Q,    &
       POSITIVE_DEFINITE_QHYD
    use mod_da_param_estimation, only: &
       DA_param_estimation_update
    implicit none

    character(len=H_LONG) :: filename
         
    integer :: i, j, k
    !---------------------------------------------------------------------------

    call PROF_rapstart('DA_Update', 1)

    !---------------------------------------------------------------------------
    ! Read observation data
    !---------------------------------------------------------------------------
    call LETKF_obs_readfile( OBS_IN_NUM, OBS_IN_FORMAT, OBS_IN_BASENAME, OBS_IN_MASKFILE )

    !---------------------------------------------------------------------------
    ! Get observation operator
    !---------------------------------------------------------------------------
    !!
    !! Compute observation operator, return the results in obsda
    !! with additional space for externally processed observations
    !!
    call ATMOS_vars_get_diagnostic( 'RH',   RH(:,:,:) )
    call ATMOS_vars_get_diagnostic( 'PREC', PREC(:,:) )

    RH(:,:,:) = RH(:,:,:) * 0.01_RP ! [%] -> [no unit]

    call LETKF_obs_operator( OBS_IN_NUM, OBS_IN_FORMAT,                            &
                             U, V, W, TEMP, PRES, QV, QC, QR, QI, QS, QG, RH, HGT, &
                             TOPO, PS, PREC, U10M, V10M, T2M, Q2M                  )

    !---------------------------------------------------------------------------
    ! Process observation data
    !---------------------------------------------------------------------------
    call LETKF_obs_initialize( OBS_IN_NUM )

    !---------------------------------------------------------------------------
    ! Output first guess
    !---------------------------------------------------------------------------
    call make_output( 'GUES' )

    !---------------------------------------------------------------------------
    ! Data Assimilation (LETKF)
    !---------------------------------------------------------------------------
    call LETKF_system( OBS_IN_NUM,              & ! [IN]
                       OBS_IN_FORMAT,           & ! [IN]
                       U   (KS:KE,IS:IE,JS:JE), & ! [INOUT]
                       V   (KS:KE,IS:IE,JS:JE), & ! [INOUT]
                       W   (KS:KE,IS:IE,JS:JE), & ! [INOUT]
                       TEMP(KS:KE,IS:IE,JS:JE), & ! [INOUT]
                       PRES(KS:KE,IS:IE,JS:JE), & ! [INOUT]
                       QV  (KS:KE,IS:IE,JS:JE), & ! [INOUT]
                       QC  (KS:KE,IS:IE,JS:JE), & ! [INOUT]
                       QR  (KS:KE,IS:IE,JS:JE), & ! [INOUT]
                       QI  (KS:KE,IS:IE,JS:JE), & ! [INOUT]
                       QS  (KS:KE,IS:IE,JS:JE), & ! [INOUT]
                       QG  (KS:KE,IS:IE,JS:JE)  ) ! [INOUT]

    !---------------------------------------------------------------------------
    ! Parameter Estimation (ETKF)
    !---------------------------------------------------------------------------
    call DA_param_estimation_update

    !---------------------------------------------------------------------------
    ! Clearing observation data structure (each step)
    !---------------------------------------------------------------------------
    call LETKF_obs_clear( OBS_IN_NUM )

    !---------------------------------------------------------------------------
    ! QV/QHYD negative filter
    !---------------------------------------------------------------------------
    if( POSITIVE_DEFINITE_Q ) then
       QV(:,:,:) = max( QV(:,:,:), 0.0_RP )
    end if
    if( POSITIVE_DEFINITE_QHYD ) then
       QC(:,:,:) = max( QC(:,:,:), 0.0_RP )
       QR(:,:,:) = max( QR(:,:,:), 0.0_RP )
       QI(:,:,:) = max( QI(:,:,:), 0.0_RP )
       QS(:,:,:) = max( QS(:,:,:), 0.0_RP )
       QG(:,:,:) = max( QG(:,:,:), 0.0_RP )
    end if

    !---------------------------------------------------------------------------
    ! overwrite the prognostic variables
    !---------------------------------------------------------------------------
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       QHYD(k,i,j,I_HC) = QC(k,i,j)
       QHYD(k,i,j,I_HR) = QR(k,i,j)
       QHYD(k,i,j,I_HI) = QI(k,i,j)
       QHYD(k,i,j,I_HS) = QS(k,i,j)
       QHYD(k,i,j,I_HG) = QG(k,i,j)
    enddo
    enddo
    enddo

    call ATMOS_PHY_MP_driver_qhyd2qtrc( KA, KS, KE, & ! [in]
                                        IA, IS, IE, & ! [in]
                                        JA, JS, JE, & ! [in]
                                        QV,         & ! [in]
                                        QHYD,       & ! [in]
                                        QTRC        ) ! [out]

    call ATMOS_THERMODYN_specific_heat( &
         KA, KS, KE, IA, IS, IE, JA, JS, JE, QA,                  & ! (in)
         QTRC(:,:,:,:),                                           & ! (in)
         TRACER_MASS(:), TRACER_R(:), TRACER_CV(:), TRACER_CP(:), & ! (in)
         Qdry(:,:,:), Rtot(:,:,:), CVtot(:,:,:), CPtot(:,:,:)     ) ! (out)

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       DENS(k,i,j) = PRES(k,i,j) / ( TEMP(k,i,j) * Rtot(k,i,j) )
       POTT(k,i,j) = TEMP(k,i,j) * ( PRE00 / PRES(k,i,j) )**(Rtot(k,i,j)/CPtot(k,i,j))
       RHOT(k,i,j) = DENS(k,i,j) * POTT(k,i,j)
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE-1
       MOMZ(k,i,j) = 0.5_RP * ( W(k,i,j)*DENS(k,i,j) + W(k+1,i,j)*DENS(k+1,i,j) )
    enddo
    enddo
    enddo
    do j = JS, JE
    do i = IS, IE
       MOMZ(KE,i,j) = W(KE,i,j) * DENS(KE,i,j)
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE-1
    do k = KS, KE
       MOMX(k,i,j) = 0.5_RP * ( U(k,i,j)*DENS(k,i,j) + U(k,i+1,j)*DENS(k,i+1,j) )
    enddo
    enddo
    enddo
    do j = JS, JE
    do k = KS, KE
       MOMX(k,IE,j) = U(k,IE,j) * DENS(k,IE,j)
    enddo
    enddo

    do j = JS, JE-1
    do i = IS, IE
    do k = KS, KE
       MOMY(k,i,j) = 0.5_RP * ( V(k,i,j)*DENS(k,i,j) + V(k,i,j+1)*DENS(k,i,j+1) )
    end do
    end do
    end do
    do i = IS, IE
    do k = KS, KE
       MOMY(k,i,JE) = V(k,i,JE) * DENS(k,i,JE)
    enddo
    enddo

    !---------------------------------------------------------------------------
    ! Output analysis
    !---------------------------------------------------------------------------
    call make_output( 'ANLS' )

    call PROF_rapend  ('DA_Update', 1)

    return
  end subroutine DA_driver_update

  subroutine make_ens_2D( &
      IMAX,    &
      JMAX,    &
      invar2d, &
      mean2d,  &
      sprd2d   )
    use mpi
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_comm_ensemble, only: &
       ENSEMBLE_world  => COMM_ENSEMBLE_world,  &
       ENSEMBLE_nprocs => COMM_ENSEMBLE_nprocs, &
       ENSEMBLE_myrank => COMM_ENSEMBLE_myrank
    use scale_statistics, only: &
       average => STATISTICS_average, &
       stddev  => STATISTICS_stddev
    implicit none

    integer,  intent(in)  :: IMAX
    integer,  intent(in)  :: JMAX
    real(RP), intent(in)  :: invar2d(:,:)
    real(RP), intent(out) :: mean2d (:,:)
    real(RP), intent(out) :: sprd2d (:,:)

    real(RP) :: send(IMAX*JMAX*ENSEMBLE_nprocs)
    real(RP) :: recv(IMAX*JMAX*ENSEMBLE_nprocs)
    real(RP) :: work(IMAX,JMAX,ENSEMBLE_nprocs)

    integer :: datasize
    integer :: i, j, n
    integer :: ierr
    !---------------------------------------------------------------------------

    datasize = IMAX * JMAX

    !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
    do j = 1, JMAX
    do i = 1, IMAX
       send( i + (j-1)*IMAX + ENSEMBLE_myrank*IMAX*JMAX ) = invar2d(i,j)
    end do
    end do

    n = ENSEMBLE_myrank * datasize + 1

    call MPI_ALLGATHER( send(n),        &
                        datasize,       &
                        datatype,       &
                        recv(1),        &
                        datasize,       &
                        datatype,       &
                        ENSEMBLE_world, &
                        ierr            )

    !$omp parallel do private(i,j,n) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
    do j = 1, JMAX
    do i = 1, IMAX
       do n = 1, ENSEMBLE_nprocs
          work(i,j,n) = recv( i + (j-1)*IMAX + (n-1)*IMAX*JMAX )
       end do
       mean2d(i,j) = average( work(i,j,:), UNDEF )
       sprd2d(i,j) = stddev ( work(i,j,:), UNDEF )
    end do
    end do

    return
  end subroutine make_ens_2D

  subroutine make_ens_3D( &
      KMAX,    &
      IMAX,    &
      JMAX,    &
      invar3d, &
      mean3d,  &
      sprd3d   )
    use mpi
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_comm_ensemble, only: &
       ENSEMBLE_world  => COMM_ENSEMBLE_world,  &
       ENSEMBLE_nprocs => COMM_ENSEMBLE_nprocs, &
       ENSEMBLE_myrank => COMM_ENSEMBLE_myrank
    use scale_statistics, only: &
       average => STATISTICS_average, &
       stddev  => STATISTICS_stddev
    implicit none

    integer,  intent(in)  :: KMAX
    integer,  intent(in)  :: IMAX
    integer,  intent(in)  :: JMAX
    real(RP), intent(in)  :: invar3d(1:KMAX,1:IMAX,1:JMAX)
    real(RP), intent(out) :: mean3d (1:KMAX,1:IMAX,1:JMAX)
    real(RP), intent(out) :: sprd3d (1:KMAX,1:IMAX,1:JMAX)

    real(RP) :: send(KMAX*IMAX*JMAX*ENSEMBLE_nprocs)
    real(RP) :: recv(KMAX*IMAX*JMAX*ENSEMBLE_nprocs)
    real(RP) :: work(KMAX,IMAX,JMAX,ENSEMBLE_nprocs)

    integer :: datasize
    integer :: k, i, j, n
    integer :: ierr
    !---------------------------------------------------------------------------

    datasize = KMAX * IMAX * JMAX

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
    do j = 1, JMAX
    do i = 1, IMAX
    do k = 1, KMAX
       send( k + (i-1)*KMAX + (j-1)*KMAX*IMAX + ENSEMBLE_myrank*KMAX*IMAX*JMAX ) = invar3d(k,i,j)
    end do
    end do
    end do

    n = ENSEMBLE_myrank * datasize + 1

    call MPI_ALLGATHER( send(n),        &
                        datasize,       &
                        datatype,       &
                        recv(1),        &
                        datasize,       &
                        datatype,       &
                        ENSEMBLE_world, &
                        ierr            )

    !$omp parallel do private(i,j,k,n) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
    do j = 1, JMAX
    do i = 1, IMAX
    do k = 1, KMAX
       do n = 1, ENSEMBLE_nprocs
          work(k,i,j,n) = recv( k + (i-1)*KMAX + (j-1)*KMAX*IMAX + (n-1)*KMAX*IMAX*JMAX )
       end do
       mean3d(k,i,j) = average( work(k,i,j,:), UNDEF )
       sprd3d(k,i,j) = stddev ( work(k,i,j,:), UNDEF )
    end do
    end do
    end do

    return
  end subroutine make_ens_3D

  subroutine make_output( PREFIX )
    use scale_atmos_grid_cartesC_index, only: &
       KMAX,       &
       IMAX,       &
       JMAX,       &
       KA, KS, KE, &
       IA, IS, IE, &
       JA, JS, JE
    use scale_file_history, only: &
       FILE_HISTORY_in
    use mod_atmos_vars, only: &
       ATMOS_vars_get_diagnostic, &
       DENS, &
       MOMX, &
       MOMY, &
       MOMZ, &
       RHOT, &
       U,    &
       V,    &
       W,    &
       TEMP, &
       PRES, &
       QV,   &
       QC,   &
       QR,   &
       QI,   &
       QS,   &
       QG
    use mod_da_vars, only: &
       DA_COMPUTE_ENS_HISTORY
    implicit none

    character(len=4), intent(in) :: PREFIX
    !---------------------------------------------------------------------------

    call ATMOS_vars_get_diagnostic( 'RH', RH(:,:,:) )

    call FILE_HISTORY_in( DENS(:,:,:), trim(PREFIX)//'_DENS', 'DENS for '//trim(PREFIX), 'kg/m3'   )
    call FILE_HISTORY_in( MOMX(:,:,:), trim(PREFIX)//'_MOMX', 'MOMX for '//trim(PREFIX), 'kg/m2/s' )
    call FILE_HISTORY_in( MOMY(:,:,:), trim(PREFIX)//'_MOMY', 'MOMY for '//trim(PREFIX), 'kg/m2/s' )
    call FILE_HISTORY_in( MOMZ(:,:,:), trim(PREFIX)//'_MOMZ', 'MOMZ for '//trim(PREFIX), 'kg/m2/s' )
    call FILE_HISTORY_in( RHOT(:,:,:), trim(PREFIX)//'_RHOT', 'RHOT for '//trim(PREFIX), 'K*kg/m3' )
    call FILE_HISTORY_in( U   (:,:,:), trim(PREFIX)//'_U',    'U for '//trim(PREFIX),    'm/s'     )
    call FILE_HISTORY_in( V   (:,:,:), trim(PREFIX)//'_V',    'V for '//trim(PREFIX),    'm/s'     )
    call FILE_HISTORY_in( W   (:,:,:), trim(PREFIX)//'_W',    'W for '//trim(PREFIX),    'm/s'     )
    call FILE_HISTORY_in( TEMP(:,:,:), trim(PREFIX)//'_T',    'TEMP for '//trim(PREFIX), 'K'       )
    call FILE_HISTORY_in( PRES(:,:,:), trim(PREFIX)//'_PRES', 'PRES for '//trim(PREFIX), 'Pa'      )
    call FILE_HISTORY_in( QV  (:,:,:), trim(PREFIX)//'_QV',   'QV for '//trim(PREFIX),   'kg/kg'   )
    call FILE_HISTORY_in( QC  (:,:,:), trim(PREFIX)//'_QC',   'QC for '//trim(PREFIX),   'kg/kg'   )
    call FILE_HISTORY_in( QR  (:,:,:), trim(PREFIX)//'_QR',   'QR for '//trim(PREFIX),   'kg/kg'   )
    call FILE_HISTORY_in( QI  (:,:,:), trim(PREFIX)//'_QI',   'QI for '//trim(PREFIX),   'kg/kg'   )
    call FILE_HISTORY_in( QS  (:,:,:), trim(PREFIX)//'_QS',   'QS for '//trim(PREFIX),   'kg/kg'   )
    call FILE_HISTORY_in( QG  (:,:,:), trim(PREFIX)//'_QG',   'QG for '//trim(PREFIX),   'kg/kg'   )
    call FILE_HISTORY_in( RH  (:,:,:), trim(PREFIX)//'_RH',   'RH for '//trim(PREFIX),   '%'       )

    !---------------------------------------------------------------------------
    ! Compute ensemble mean/spread
    !---------------------------------------------------------------------------
    if( DA_COMPUTE_ENS_HISTORY ) then
      call make_ens_3D( KMAX, IMAX, JMAX, DENS(KS:KE,IS:IE,JS:JE), ENS_mean_DENS(KS:KE,IS:IE,JS:JE), ENS_sprd_DENS(KS:KE,IS:IE,JS:JE) )
      call make_ens_3D( KMAX, IMAX, JMAX, MOMX(KS:KE,IS:IE,JS:JE), ENS_mean_MOMX(KS:KE,IS:IE,JS:JE), ENS_sprd_MOMX(KS:KE,IS:IE,JS:JE) )
      call make_ens_3D( KMAX, IMAX, JMAX, MOMY(KS:KE,IS:IE,JS:JE), ENS_mean_MOMY(KS:KE,IS:IE,JS:JE), ENS_sprd_MOMY(KS:KE,IS:IE,JS:JE) )
      call make_ens_3D( KMAX, IMAX, JMAX, MOMZ(KS:KE,IS:IE,JS:JE), ENS_mean_MOMZ(KS:KE,IS:IE,JS:JE), ENS_sprd_MOMZ(KS:KE,IS:IE,JS:JE) )
      call make_ens_3D( KMAX, IMAX, JMAX, RHOT(KS:KE,IS:IE,JS:JE), ENS_mean_RHOT(KS:KE,IS:IE,JS:JE), ENS_sprd_RHOT(KS:KE,IS:IE,JS:JE) )
      call make_ens_3D( KMAX, IMAX, JMAX, U   (KS:KE,IS:IE,JS:JE), ENS_mean_U   (KS:KE,IS:IE,JS:JE), ENS_sprd_U   (KS:KE,IS:IE,JS:JE) )
      call make_ens_3D( KMAX, IMAX, JMAX, V   (KS:KE,IS:IE,JS:JE), ENS_mean_V   (KS:KE,IS:IE,JS:JE), ENS_sprd_V   (KS:KE,IS:IE,JS:JE) )
      call make_ens_3D( KMAX, IMAX, JMAX, W   (KS:KE,IS:IE,JS:JE), ENS_mean_W   (KS:KE,IS:IE,JS:JE), ENS_sprd_W   (KS:KE,IS:IE,JS:JE) )
      call make_ens_3D( KMAX, IMAX, JMAX, TEMP(KS:KE,IS:IE,JS:JE), ENS_mean_TEMP(KS:KE,IS:IE,JS:JE), ENS_sprd_TEMP(KS:KE,IS:IE,JS:JE) )
      call make_ens_3D( KMAX, IMAX, JMAX, PRES(KS:KE,IS:IE,JS:JE), ENS_mean_PRES(KS:KE,IS:IE,JS:JE), ENS_sprd_PRES(KS:KE,IS:IE,JS:JE) )
      call make_ens_3D( KMAX, IMAX, JMAX, QV  (KS:KE,IS:IE,JS:JE), ENS_mean_QV  (KS:KE,IS:IE,JS:JE), ENS_sprd_QV  (KS:KE,IS:IE,JS:JE) )
      call make_ens_3D( KMAX, IMAX, JMAX, QC  (KS:KE,IS:IE,JS:JE), ENS_mean_QC  (KS:KE,IS:IE,JS:JE), ENS_sprd_QC  (KS:KE,IS:IE,JS:JE) )
      call make_ens_3D( KMAX, IMAX, JMAX, QR  (KS:KE,IS:IE,JS:JE), ENS_mean_QR  (KS:KE,IS:IE,JS:JE), ENS_sprd_QR  (KS:KE,IS:IE,JS:JE) )
      call make_ens_3D( KMAX, IMAX, JMAX, QI  (KS:KE,IS:IE,JS:JE), ENS_mean_QI  (KS:KE,IS:IE,JS:JE), ENS_sprd_QI  (KS:KE,IS:IE,JS:JE) )
      call make_ens_3D( KMAX, IMAX, JMAX, QS  (KS:KE,IS:IE,JS:JE), ENS_mean_QS  (KS:KE,IS:IE,JS:JE), ENS_sprd_QS  (KS:KE,IS:IE,JS:JE) )
      call make_ens_3D( KMAX, IMAX, JMAX, QG  (KS:KE,IS:IE,JS:JE), ENS_mean_QG  (KS:KE,IS:IE,JS:JE), ENS_sprd_QG  (KS:KE,IS:IE,JS:JE) )
      call make_ens_3D( KMAX, IMAX, JMAX, RH  (KS:KE,IS:IE,JS:JE), ENS_mean_RH  (KS:KE,IS:IE,JS:JE), ENS_sprd_RH  (KS:KE,IS:IE,JS:JE) )

      call FILE_HISTORY_in( ENS_mean_DENS(:,:,:), trim(PREFIX)//'_mean_DENS', 'Ensemble mean of DENS for '//trim(PREFIX), 'kg/m3'   )
      call FILE_HISTORY_in( ENS_mean_MOMX(:,:,:), trim(PREFIX)//'_mean_MOMX', 'Ensemble mean of MOMX for '//trim(PREFIX), 'kg/m2/s' )
      call FILE_HISTORY_in( ENS_mean_MOMY(:,:,:), trim(PREFIX)//'_mean_MOMY', 'Ensemble mean of MOMY for '//trim(PREFIX), 'kg/m2/s' )
      call FILE_HISTORY_in( ENS_mean_MOMZ(:,:,:), trim(PREFIX)//'_mean_MOMZ', 'Ensemble mean of MOMZ for '//trim(PREFIX), 'kg/m2/s' )
      call FILE_HISTORY_in( ENS_mean_RHOT(:,:,:), trim(PREFIX)//'_mean_RHOT', 'Ensemble mean of RHOT for '//trim(PREFIX), 'K*kg/m3' )
      call FILE_HISTORY_in( ENS_mean_U   (:,:,:), trim(PREFIX)//'_mean_U',    'Ensemble mean of U for '//trim(PREFIX),    'm/s'     )
      call FILE_HISTORY_in( ENS_mean_V   (:,:,:), trim(PREFIX)//'_mean_V',    'Ensemble mean of V for '//trim(PREFIX),    'm/s'     )
      call FILE_HISTORY_in( ENS_mean_W   (:,:,:), trim(PREFIX)//'_mean_W',    'Ensemble mean of W for '//trim(PREFIX),    'm/s'     )
      call FILE_HISTORY_in( ENS_mean_TEMP(:,:,:), trim(PREFIX)//'_mean_T',    'Ensemble mean of TEMP for '//trim(PREFIX), 'K'       )
      call FILE_HISTORY_in( ENS_mean_PRES(:,:,:), trim(PREFIX)//'_mean_PRES', 'Ensemble mean of PRES for '//trim(PREFIX), 'Pa'      )
      call FILE_HISTORY_in( ENS_mean_QV  (:,:,:), trim(PREFIX)//'_mean_QV',   'Ensemble mean of QV for '//trim(PREFIX),   'kg/kg'   )
      call FILE_HISTORY_in( ENS_mean_QC  (:,:,:), trim(PREFIX)//'_mean_QC',   'Ensemble mean of QC for '//trim(PREFIX),   'kg/kg'   )
      call FILE_HISTORY_in( ENS_mean_QR  (:,:,:), trim(PREFIX)//'_mean_QR',   'Ensemble mean of QR for '//trim(PREFIX),   'kg/kg'   )
      call FILE_HISTORY_in( ENS_mean_QI  (:,:,:), trim(PREFIX)//'_mean_QI',   'Ensemble mean of QI for '//trim(PREFIX),   'kg/kg'   )
      call FILE_HISTORY_in( ENS_mean_QS  (:,:,:), trim(PREFIX)//'_mean_QS',   'Ensemble mean of QS for '//trim(PREFIX),   'kg/kg'   )
      call FILE_HISTORY_in( ENS_mean_QG  (:,:,:), trim(PREFIX)//'_mean_QG',   'Ensemble mean of QG for '//trim(PREFIX),   'kg/kg'   )
      call FILE_HISTORY_in( ENS_mean_RH  (:,:,:), trim(PREFIX)//'_mean_RH',   'Ensemble mean of RH for '//trim(PREFIX),   '%'       )

      call FILE_HISTORY_in( ENS_sprd_DENS(:,:,:), trim(PREFIX)//'_sprd_DENS', 'Ensemble spread of DENS for '//trim(PREFIX), 'kg/m3'   )
      call FILE_HISTORY_in( ENS_sprd_MOMX(:,:,:), trim(PREFIX)//'_sprd_MOMX', 'Ensemble spread of MOMX for '//trim(PREFIX), 'kg/m2/s' )
      call FILE_HISTORY_in( ENS_sprd_MOMY(:,:,:), trim(PREFIX)//'_sprd_MOMY', 'Ensemble spread of MOMY for '//trim(PREFIX), 'kg/m2/s' )
      call FILE_HISTORY_in( ENS_sprd_MOMZ(:,:,:), trim(PREFIX)//'_sprd_MOMZ', 'Ensemble spread of MOMZ for '//trim(PREFIX), 'kg/m2/s' )
      call FILE_HISTORY_in( ENS_sprd_RHOT(:,:,:), trim(PREFIX)//'_sprd_RHOT', 'Ensemble spread of RHOT for '//trim(PREFIX), 'K*kg/m3' )
      call FILE_HISTORY_in( ENS_sprd_U   (:,:,:), trim(PREFIX)//'_sprd_U',    'Ensemble spread of U for '//trim(PREFIX),    'm/s'     )
      call FILE_HISTORY_in( ENS_sprd_V   (:,:,:), trim(PREFIX)//'_sprd_V',    'Ensemble spread of V for '//trim(PREFIX),    'm/s'     )
      call FILE_HISTORY_in( ENS_sprd_W   (:,:,:), trim(PREFIX)//'_sprd_W',    'Ensemble spread of W for '//trim(PREFIX),    'm/s'     )
      call FILE_HISTORY_in( ENS_sprd_TEMP(:,:,:), trim(PREFIX)//'_sprd_T',    'Ensemble spread of TEMP for '//trim(PREFIX), 'K'       )
      call FILE_HISTORY_in( ENS_sprd_PRES(:,:,:), trim(PREFIX)//'_sprd_PRES', 'Ensemble spread of PRES for '//trim(PREFIX), 'Pa'      )
      call FILE_HISTORY_in( ENS_sprd_QV  (:,:,:), trim(PREFIX)//'_sprd_QV',   'Ensemble spread of QV for '//trim(PREFIX),   'kg/kg'   )
      call FILE_HISTORY_in( ENS_sprd_QC  (:,:,:), trim(PREFIX)//'_sprd_QC',   'Ensemble spread of QC for '//trim(PREFIX),   'kg/kg'   )
      call FILE_HISTORY_in( ENS_sprd_QR  (:,:,:), trim(PREFIX)//'_sprd_QR',   'Ensemble spread of QR for '//trim(PREFIX),   'kg/kg'   )
      call FILE_HISTORY_in( ENS_sprd_QI  (:,:,:), trim(PREFIX)//'_sprd_QI',   'Ensemble spread of QI for '//trim(PREFIX),   'kg/kg'   )
      call FILE_HISTORY_in( ENS_sprd_QS  (:,:,:), trim(PREFIX)//'_sprd_QS',   'Ensemble spread of QS for '//trim(PREFIX),   'kg/kg'   )
      call FILE_HISTORY_in( ENS_sprd_QG  (:,:,:), trim(PREFIX)//'_sprd_QG',   'Ensemble spread of QG for '//trim(PREFIX),   'kg/kg'   )
      call FILE_HISTORY_in( ENS_sprd_RH  (:,:,:), trim(PREFIX)//'_sprd_RH',   'Ensemble spread of RH for '//trim(PREFIX),   '%'       )
    end if

    return
  end subroutine make_output

end module mod_da_driver
