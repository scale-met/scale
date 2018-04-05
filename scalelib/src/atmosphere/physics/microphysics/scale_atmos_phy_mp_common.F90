!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics - Common
!!
!! @par Description
!!          Common module for Cloud Microphysics
!!          Sedimentation/Precipitation and Saturation adjustment
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2012-12-23 (H.Yashiro)  [new]
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_phy_mp_common
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_atmos_grid_cartesC_index
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_MP_negative_fixer
  public :: ATMOS_PHY_MP_saturation_adjustment
  public :: ATMOS_PHY_MP_precipitation
  public :: ATMOS_PHY_MP_precipitation_momentum

  interface ATMOS_PHY_MP_saturation_adjustment
     module procedure ATMOS_PHY_MP_saturation_adjustment_3D
  end interface ATMOS_PHY_MP_saturation_adjustment

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> ATMOS_PHY_MP_negative_fixer
  !! negative fixer
  !<
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_negative_fixer( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, QHA, &
       limit_negative, &
       DENS, QV, QTRC )
    use scale_prc, only: &
       PRC_myrank, &
       PRC_abort
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    integer, intent(in) :: QHA

    real(RP), intent(in)    :: limit_negative

    real(RP), intent(inout) :: DENS(KA,IA,JA)
    real(RP), intent(inout) :: QV  (KA,IA,JA)
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QHA)

    real(RP) :: diffq
    real(RP) :: diffq_check(KA,IA,JA)
    real(RP) :: diffq_min

    integer  :: k, i, j, iq
    !---------------------------------------------------------------------------

    call PROF_rapstart('MP_filter', 3)

    !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
    !$omp private(i,j,k, &
    !$omp         iq,diffq) &
    !$omp shared(KS,KE,IS,IE,JS,JE,QHA, &
    !$omp        QV,QTRC,DENS,diffq_check)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE

       diffq = 0.0_RP
       do iq = 1, QHA
          ! total hydrometeor (before correction)
          diffq = diffq + QTRC(k,i,j,iq)
          ! remove negative value of hydrometeors (mass)
          QTRC(k,i,j,iq) = max( QTRC(k,i,j,iq), 0.0_RP )
       enddo

       do iq = 1, QHA
          ! difference between before and after correction
          diffq = diffq - QTRC(k,i,j,iq)
       enddo

       ! Compensate for the lack of hydrometeors by the water vapor
       QV(k,i,j) = QV(k,i,j) + diffq
       diffq_check(k,i,j) = diffq

       ! TODO: We have to consider energy conservation (but very small)

       ! remove negative value of water vapor (mass)
       diffq = QV(k,i,j)
       QV(k,i,j) = max( QV(k,i,j), 0.0_RP )
       diffq = diffq - QV(k,i,j)

       ! Apply correction to total density
       ! TODO: We have to consider energy conservation (but very small)
       DENS(k,i,j) = DENS(k,i,j) * ( 1.0_RP - diffq ) ! diffq is negative

    enddo
    enddo
    enddo

    diffq_min = minval( diffq_check(KS:KE,ISB:IEB,JSB:JEB) )

    if (       abs(limit_negative) > 0.0_RP         &
         .AND. abs(limit_negative) < abs(diffq_min) ) then
       LOG_ERROR("ATMOS_PHY_MP_negative_fixer",*) 'large negative is found. rank = ', PRC_myrank

       do j = JS, JE
       do i = IS, IE
       do k = KS,  KE
          if (     abs(limit_negative) < abs(diffq_check(k,i,j)) &
              .OR. abs(QV(k,i,j)     ) < abs(diffq_check(k,i,j)) ) then
             LOG_ERROR_CONT(*) 'k,i,j,value(QHYD,QV) = ', k, i, j, diffq_check(k,i,j), QV(k,i,j)
          endif
       enddo
       enddo
       enddo
       LOG_ERROR_CONT(*) 'total negative hydrometeor < ', abs(limit_negative)

       call PRC_abort
    endif

    call PROF_rapend('MP_filter', 3)

    return
  end subroutine ATMOS_PHY_MP_negative_fixer

  subroutine ATMOS_PHY_MP_negative_fixer_obsolute( &
       DENS,          &
       RHOT,          &
       QTRC,          &
       I_QV,          &
       limit_negative )
    use scale_prc, only: &
       PRC_myrank, &
       PRC_abort
    use scale_atmos_hydrometeor, only: &
       QHS, &
       QHE
    implicit none

    real(RP), intent(inout) :: DENS(KA,IA,JA)
    real(RP), intent(inout) :: RHOT(KA,IA,JA)
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA)
    integer,  intent(in)    :: I_QV
    real(RP), intent(in)    :: limit_negative

    real(RP) :: diffq
    real(RP) :: diffq_check(KA,IA,JA)
    real(RP) :: diffq_min

    integer  :: k, i, j, iq
    !---------------------------------------------------------------------------

    call PROF_rapstart('MP_filter', 3)

    !$omp parallel do default(none) private(i,j,k,iq,diffq) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(JSB,JEB,ISB,IEB,KS,KE,QHS,QHE,QTRC,DENS,RHOT,I_QV,diffq_check)
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS,  KE

       diffq = 0.0_RP
       do iq = QHS, QHE
          ! total hydrometeor (before correction)
          diffq = diffq + QTRC(k,i,j,iq)
          ! remove negative value of hydrometeors (mass)
          QTRC(k,i,j,iq) = max( QTRC(k,i,j,iq), 0.0_RP )
       enddo

       do iq = QHS, QHE
          ! difference between before and after correction
          diffq = diffq - QTRC(k,i,j,iq)
       enddo

       ! Compensate for the lack of hydrometeors by the water vapor
       QTRC(k,i,j,I_QV) = QTRC(k,i,j,I_QV) + diffq
       diffq_check(k,i,j) = diffq

       ! TODO: We have to consider energy conservation (but very small)

       ! remove negative value of water vapor (mass)
       diffq = QTRC(k,i,j,I_QV)
       QTRC(k,i,j,I_QV) = max( QTRC(k,i,j,I_QV), 0.0_RP )
       diffq = diffq - QTRC(k,i,j,I_QV)

       ! Apply correction to total density
       ! TODO: We have to consider energy conservation (but very small)
       DENS(k,i,j) = DENS(k,i,j) * ( 1.0_RP - diffq ) ! diffq is negative
       RHOT(k,i,j) = RHOT(k,i,j) * ( 1.0_RP - diffq )

    enddo
    enddo
    enddo

    diffq_min = minval( diffq_check(KS:KE,ISB:IEB,JSB:JEB) )

    if (       abs(limit_negative) > 0.0_RP         &
         .AND. abs(limit_negative) < abs(diffq_min) ) then
       LOG_ERROR("ATMOS_PHY_MP_negative_fixer_obsolute",*) 'large negative is found. rank = ', PRC_myrank

       do j = JSB, JEB
       do i = ISB, IEB
       do k = KS,  KE
          if (     abs(limit_negative)   < abs(diffq_check(k,i,j)) &
              .OR. abs(QTRC(k,i,j,I_QV)) < abs(diffq_check(k,i,j)) ) then
             LOG_ERROR_CONT(*) 'k,i,j,value(QHYD,QV) = ', k, i, j, diffq_check(k,i,j), QTRC(k,i,j,I_QV)
          endif
       enddo
       enddo
       enddo
       LOG_ERROR_CONT(*) 'total negative hydrometeor < ', abs(limit_negative)

       call PRC_abort
    endif

    call PROF_rapend('MP_filter', 3)

    return
  end subroutine ATMOS_PHY_MP_negative_fixer_obsolute

  !-----------------------------------------------------------------------------
  !> Saturation adjustment
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_saturation_adjustment_3D( &
         KA, KS, KE, IA, IS, IE, JA, JS, JE, &
         DENS,         &
         flag_liquid,  &
         TEMP,         &
         QV, QC, QI,   &
         CPtot, CVtot, &
         RHOE_d        )
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_hydrometeor, only: &
       CP_VAPOR, &
       CP_WATER, &
       CP_ICE,   &
       CV_VAPOR, &
       CV_WATER, &
       CV_ICE,   &
       LHV,   &
       LHF
    use scale_atmos_saturation, only: &
       ATMOS_SATURATION_moist_conversion_dens_all, &
       ATMOS_SATURATION_moist_conversion_dens_liq
    implicit none

    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in) :: DENS(KA,IA,JA)
    logical , intent(in) :: flag_liquid !> use scheme only for the liquid water?

    real(RP), intent(inout) :: TEMP (KA,IA,JA)
    real(RP), intent(inout) :: QV   (KA,IA,JA)
    real(RP), intent(inout) :: QC   (KA,IA,JA)
    real(RP), intent(inout) :: QI   (KA,IA,JA)
    real(RP), intent(inout) :: CPtot(KA,IA,JA)
    real(RP), intent(inout) :: CVtot(KA,IA,JA)

    real(RP), intent(out) :: RHOE_d(KA,IA,JA)

    ! working
    real(RP) :: QV1
    real(RP) :: QC1
    real(RP) :: QI1

    real(RP) :: Emoist ! moist internal energy

    logical :: converged, error

    integer :: k, i, j
    !---------------------------------------------------------------------------

    call PROF_rapstart('MP_Saturation_adjustment', 2)

    error = .false.

    if ( flag_liquid ) then ! warm rain

       !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
       !$omp shared(KS,KE,IS,IE,JS,JE, &
       !$omp        CP_VAPOR,CP_WATER,CV_VAPOR,CV_WATER,LHV,LHF, &
       !$omp        DENS,QV,QC,TEMP,CPtot,CVtot,RHOE_d,error) &
       !$omp private(i,j,k, &
       !$omp         QV1,QC1,Emoist,converged)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE

          QV1 = QV(k,i,j)
          QC1 = QC(k,i,j)

          Emoist = TEMP(k,i,j) * CVtot(k,i,j) &
                 + QV1 * LHV

          call ATMOS_SATURATION_moist_conversion_dens_liq( DENS(k,i,j), Emoist,        & ! [IN]
                                                           TEMP(k,i,j), QV1, QC1,      & ! [INOUT]
                                                           CPtot(k,i,j), CVtot(k,i,j), & ! [INOUT]
                                                           converged                   ) ! [OUT]

          if ( .NOT. converged ) then
             LOG_ERROR("ATMOS_PHY_MP_saturation_adjustment_3D",*) 'moist_conversion not converged! ', k,i,j
             error = .true.
             exit
          endif

          RHOE_d(k,i,j) = - LHV * ( QV1 - QV(k,i,j) ) * DENS(k,i,j)

          QV(k,i,j) = QV1
          QC(k,i,j) = QC1

       end do
       end do
       end do

       if ( error ) call PRC_abort

    else ! cold rain

       !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
       !$omp shared (KS,KE,IS,IE,JS,JE, &
       !$omp         CP_VAPOR,CP_WATER,CP_ICE,CV_VAPOR,CV_WATER,CV_ICE,LHV,LHF, &
       !$omp         DENS,QV,QC,QI,TEMP,CPtot,CVtot,RHOE_d,error) &
       !$omp private(i,j,k, &
       !$omp         QV1,QC1,QI1,Emoist,converged)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          QV1 = QV(k,i,j)
          QC1 = QC(k,i,j)
          QI1 = QI(k,i,j)

          Emoist = TEMP(k,i,j) * CVtot(k,i,j) &
                 + QV1 * LHV &
                 - QI1 * LHF

          call ATMOS_SATURATION_moist_conversion_dens_all( DENS(k,i,j), Emoist,        & ! [IN]
                                                           TEMP(k,i,j), QV1, QC1, QI1, & ! [INOUT]
                                                           CPtot(k,i,j), CVtot(k,i,j), & ! [INOUT]
                                                           converged                   ) ! [OUT]

          if ( .NOT. converged ) then
             LOG_ERROR("ATMOS_PHY_MP_saturation_adjustment_3D",*) 'moist_conversion not converged! ', k,i,j
             error = .true.
             exit
          endif

          RHOE_d(k,i,j) = ( - LHV * ( QV1 - QV(k,i,j) ) &
                            + LHF * ( QI1 - QI(k,i,j) ) ) * DENS(k,i,j)

          QV(k,i,j) = QV1
          QC(k,i,j) = QC1
          QI(k,i,j) = QI1

       end do
       end do
       end do

       if ( error ) call PRC_abort
    endif

    call PROF_rapend  ('MP_Saturation_adjustment', 2)

    return
  end subroutine ATMOS_PHY_MP_saturation_adjustment_3D

  !-----------------------------------------------------------------------------
  !> ATMOS_PHY_MP_precipitation
  !! precipitation transport
  !<
!OCL SERIAL
  subroutine ATMOS_PHY_MP_precipitation( &
       KA, KS, KE, QHA, QLA, QIA, &
       TEMP, vterm, FDZ, RCDZ, dt,     &
       i, j,                           &
       DENS, RHOQ, CPtot, CVtot, RHOE, &
       mflx, sflx                      )
    use scale_const, only: &
       GRAV  => CONST_GRAV
    use scale_atmos_hydrometeor, only: &
       CP_WATER, &
       CP_ICE, &
       CV_WATER, &
       CV_ICE
    implicit none
    integer,  intent(in) :: KA, KS, KE
    integer,  intent(in) :: QHA, QLA, QIA

    real(RP), intent(in) :: TEMP (KA)
    real(RP), intent(in) :: vterm(KA,QHA) ! terminal velocity of cloud mass
    real(RP), intent(in) :: FDZ  (KA)
    real(RP), intent(in) :: RCDZ (KA)
    real(DP), intent(in) :: dt
    integer,  intent(in) :: i, j         ! for debug

    real(RP), intent(inout) :: DENS (KA)
    real(RP), intent(inout) :: RHOQ (KA,QHA)
    real(RP), intent(inout) :: CPtot(KA)
    real(RP), intent(inout) :: CVtot(KA)
    real(RP), intent(inout) :: RHOE (KA)

    real(RP), intent(out)   :: mflx (KA)
    real(RP), intent(out)   :: sflx (2) !> 1: rain, 2: snow

    real(RP) :: qflx(KA)
    real(RP) :: eflx(KA)
    real(RP) :: RHOCP(KA)
    real(RP) :: RHOCV(KA)
    real(RP) :: dDENS
    real(RP) :: CP, CV

    integer  :: k, iq
    !---------------------------------------------------------------------------

    ! tracer/energy transport by falldown
    ! 1st order upwind, forward euler, velocity is always negative

    mflx(:) = 0.0_RP
    sflx(:) = 0.0_RP
    qflx(KE) = 0.0_RP
    eflx(KE) = 0.0_RP

    do k = KS, KE
       RHOCP(k) = CPtot(k) * DENS(k)
       RHOCV(k) = CVtot(k) * DENS(k)
    end do

    do iq = 1, QHA

       !--- mass flux for each tracer, upwind with vel < 0
       do k = KS-1, KE-1
          qflx(k) = vterm(k+1,iq) * RHOQ(k+1,iq)
       enddo

       !--- update falling tracer
       do k  = KS, KE
          rhoq(k,iq) = rhoq(k,iq) - dt * ( qflx(k) - qflx(k-1) ) * RCDZ(k)
       enddo ! falling (water mass & number) tracer

       ! QTRC(iq; iq>QLA+QLI) is not mass tracer, such as number density
       if ( iq > QLA + QIA ) cycle

       do k = KS-1, KE-1
          mflx(k) = mflx(k) + qflx(k)
       end do

       if ( iq > QLA ) then ! ice water
          CP = CP_ICE
          CV = CV_ICE
          sflx(2) = sflx(2) + qflx(KS-1)
       else                 ! liquid water
          CP = CP_WATER
          CV = CV_WATER
          sflx(1) = sflx(1) + qflx(KS-1)
       end if

       !--- update density
       do k = KS, KE
          dDENS = - ( qflx(k) - qflx(k-1) ) * RCDZ(k) * dt
          RHOCP(k) = RHOCP(k) + CP * dDENS
          RHOCV(k) = RHOCV(k) + CV * dDENS
          DENS(k) = DENS(k) + dDENS
       end do

       ! internal energy flux
       do k = KS-1, KE-1
          eflx(k) = qflx(k) * TEMP(k+1) * CV &
                  + qflx(k) * FDZ(k) * GRAV               ! potential energy
       end do
       !--- update internal energy
       do k = KS, KE
          RHOE(k) = RHOE(k) - ( eflx(k) - eflx(k-1) ) * RCDZ(k) * dt
       end do

    end do

    do k = KS, KE
       CPtot(k) = RHOCP(k) / DENS(k)
       CVtot(k) = RHOCV(k) / DENS(k)
    end do

    return
  end subroutine ATMOS_PHY_MP_precipitation
  !-----------------------------------------------------------------------------
  !> ATMOS_PHY_MP_precipitation_transfer
  !! precipitation transport
  !<
!OCL SERIAL
  subroutine ATMOS_PHY_MP_precipitation_momentum( &
       KA, KS, KE, &
       DENS, MOMZ, U, V, mflx, &
       RCDZ, RFDZ,             &
       MOMZ_t, RHOU_t, RHOV_t  )
    implicit none
    integer,  intent(in) :: KA, KS, KE

    real(RP), intent(in) :: DENS(KA)
    real(RP), intent(in) :: MOMZ(KA)
    real(RP), intent(in) :: U   (KA)
    real(RP), intent(in) :: V   (KA)
    real(RP), intent(in) :: mflx(KA)
    real(RP), intent(in) :: RCDZ(KA)
    real(RP), intent(in) :: RFDZ(KA)

    real(RP), intent(out) :: MOMZ_t(KA)
    real(RP), intent(out) :: RHOU_t(KA)
    real(RP), intent(out) :: RHOV_t(KA)

    real(RP) :: flx(KA)

    integer  :: k
    integer  :: iq, iqa
    !---------------------------------------------------------------------------

    flx(KE) = 0.0_RP

    !--- momentum z (half level)
    do k = KS, KE-2
       flx(k) = ( mflx(k) + mflx(k-1) ) * MOMZ(k+1) / ( DENS(k+2) + DENS(k+1) )
    enddo
    flx(KE-1) = 0.0_RP
    do k  = KS, KE-1
       MOMZ_t(k) = - ( flx(k+1) - flx(k) ) * RFDZ(k)
    enddo
    MOMZ_t(KE) = 0.0_RP

    !--- momentum x
    do k = KS-1, KE-1
       flx(k) = mflx(k) * U(k+1)
    enddo
    do k = KS, KE
       RHOU_t(k) = - ( flx(k) - flx(k-1) ) * RCDZ(k)
    enddo

    !--- momentum y
    do k = KS-1, KE-1
       flx(k) = mflx(k) * V(k+1)
    enddo
    do k = KS, KE
       RHOV_t(k) = - ( flx(k) - flx(k-1) ) * RCDZ(k)
    enddo

    return
  end subroutine ATMOS_PHY_MP_precipitation_momentum

  subroutine ATMOS_PHY_MP_precipitation_obsolute( &
       QA_MP,   &
       QS_MP,   &
       qflx,    &
       vterm,   &
       DENS,    &
       MOMZ,    &
       MOMX,    &
       MOMY,    &
       RHOE,    &
       QTRC,    &
       temp,    &
       CVq,     &
       dt,      &
       vt_fixed )
    use scale_const, only: &
       GRAV  => CONST_GRAV
    use scale_atmos_grid_cartesC_real, only: &
       REAL_CZ => ATMOS_GRID_CARTESC_REAL_CZ, &
       REAL_FZ => ATMOS_GRID_CARTESC_REAL_FZ
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    integer,  intent(in)    :: QA_MP
    integer,  intent(in)    :: QS_MP
    real(RP), intent(out)   :: qflx (KA,IA,JA,QA_MP-1)
    real(RP), intent(inout) :: vterm(KA,IA,JA,QA_MP-1) ! terminal velocity of cloud mass
    real(RP), intent(inout) :: DENS (KA,IA,JA)
    real(RP), intent(inout) :: MOMZ (KA,IA,JA)
    real(RP), intent(inout) :: MOMX (KA,IA,JA)
    real(RP), intent(inout) :: MOMY (KA,IA,JA)
    real(RP), intent(inout) :: RHOE (KA,IA,JA)
    real(RP), intent(inout) :: QTRC (KA,IA,JA,QA)
    real(RP), intent(in)    :: temp (KA,IA,JA)
    real(RP), intent(in)    :: CVq  (QA)
    real(DP), intent(in)    :: dt
    logical,  intent(in), optional :: vt_fixed

    real(RP) :: rhoq  (KA,IA,JA,QA) ! rho * q before precipitation
    real(RP) :: eflx  (KA,IA,JA)
    real(RP) :: rfdz  (KA,IA,JA)
    real(RP) :: rcdz  (KA,IA,JA)
    real(RP) :: rcdz_u(KA,IA,JA)
    real(RP) :: rcdz_v(KA,IA,JA)

    integer  :: k, i, j
    integer  :: iq, iqa
    logical  :: vt_fixed_
    !---------------------------------------------------------------------------

    call PROF_rapstart('MP_Precipitation', 2)

    if ( present(vt_fixed) ) then
       vt_fixed_ = vt_fixed
    else
       vt_fixed_ = .false.
    end if

    do iq = 1, QA_MP-1
       iqa = QS_MP + iq

       if( TRACER_MASS(iqa) == 0.0_RP ) cycle

       if( .NOT. vt_fixed_ ) call COMM_vars8( vterm(:,:,:,iq), iq )

       call COMM_vars8( QTRC(:,:,:,iqa), QA_MP+iq )
    enddo

    ! tracer/energy transport by falldown
    ! 1st order upwind, forward euler, velocity is always negative

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
       rfdz(KS-1,i,j) = 1.0_RP / ( REAL_CZ(KS,i,j) - REAL_FZ(KS-1,i,j) )
       do k = KS, KE
          rfdz  (k,i,j) = 1.0_RP / ( REAL_CZ(k+1,i,j) - REAL_CZ(k  ,i,j) )
          rcdz  (k,i,j) = 1.0_RP / ( REAL_FZ(k  ,i,j) - REAL_FZ(k-1,i,j) )

          rcdz_u(k,i,j) = 2.0_RP / ( ( REAL_FZ(k,i+1,j) - REAL_FZ(k-1,i+1,j) ) &
                                   + ( REAL_FZ(k,i  ,j) - REAL_FZ(k-1,i  ,j) ) )
          rcdz_v(k,i,j) = 2.0_RP / ( ( REAL_FZ(k,i,j+1) - REAL_FZ(k-1,i,j+1) ) &
                                   + ( REAL_FZ(k,i,j  ) - REAL_FZ(k-1,i,j  ) ) )
       enddo
    enddo
    enddo

    do iq = 1, QA_MP-1
       iqa = QS_MP + iq
       if ( TRACER_MASS(iqa) == 0.0_RP ) cycle

       if ( .not. vt_fixed_ ) then
          call COMM_wait( vterm(:,:,:,iq), iq )
       endif
       call COMM_wait( QTRC(:,:,:,iqa), QA_MP+iq )

       !$omp parallel do default(none)                                                   &
       !$omp shared(JS,JE,IS,IE,KS,KE,qflx,iq,vterm,DENS,QTRC,iqa,eflx,temp,CVq,RHOE,dt) &
       !$omp shared(rcdz,GRAV,rfdz,MOMZ,MOMX,rcdz_u,MOMY,rcdz_v)                         &
       !$omp private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j  = JS, JE
       do i  = IS, IE

          !--- mass flux for each mass tracer, upwind with vel < 0
          do k  = KS-1, KE-1
             qflx(k,i,j,iq) = vterm(k+1,i,j,iq) * DENS(k+1,i,j) * QTRC(k+1,i,j,iqa)
          enddo
          qflx(KE,i,j,iq) = 0.0_RP

          !--- internal energy
          eflx(KS-1,i,j) = qflx(KS-1,i,j,iq) * temp(KS,i,j) * CVq(iqa)
          do k  = KS, KE-1
             eflx(k,i,j) = qflx(k,i,j,iq) * temp(k+1,i,j) * CVq(iqa)
             RHOE(k,i,j) = RHOE(k,i,j) - dt * ( eflx(k,i,j) - eflx(k-1,i,j) ) * rcdz(k,i,j)
          enddo
          eflx(KE,i,j) = 0.0_RP
          RHOE(KE,i,j) = RHOE(KE,i,j) - dt * ( - eflx(KE-1,i,j) ) * rcdz(KE,i,j)

          !--- potential energy
          eflx(KS-1,i,j) = qflx(KS-1,i,j,iq) * GRAV / rfdz(KS-1,i,j)
          do k  = KS, KE-1
             eflx(k,i,j) = qflx(k,i,j,iq) * GRAV / rfdz(k,i,j)
             RHOE(k,i,j) = RHOE(k,i,j) - dt * ( eflx(k,i,j) - eflx(k-1,i,j) ) * rcdz(k,i,j)
          enddo
          RHOE(KE,i,j) = RHOE(KE,i,j) - dt * ( - eflx(KE-1,i,j) ) * rcdz(KE,i,j)

          !--- momentum z (half level)
          do k  = KS-1, KE-1
             eflx(k,i,j) = 0.25_RP * ( vterm(k+1,i,j,iq ) + vterm(k,i,j,iq ) ) &
                                   * ( QTRC (k+1,i,j,iqa) + QTRC (k,i,j,iqa) ) &
                                   * MOMZ(k,i,j)
          enddo
          do k  = KS, KE-1
             MOMZ(k,i,j) = MOMZ(k,i,j) - dt * ( eflx(k+1,i,j) - eflx(k,i,j) ) * rfdz(k,i,j)
          enddo

          !--- momentum x
          eflx(KS-1,i,j) = 0.25_RP * ( vterm(KS,i,j,iq ) + vterm(KS,i+1,j,iq ) ) &
                                   * ( QTRC (KS,i,j,iqa) + QTRC (KS,i+1,j,iqa) ) &
                                   * MOMX(KS,i,j)
          do k  = KS, KE-1
             eflx(k,i,j) = 0.25_RP * ( vterm(k+1,i,j,iq ) + vterm(k+1,i+1,j,iq ) ) &
                                   * ( QTRC (k+1,i,j,iqa) + QTRC (k+1,i+1,j,iqa) ) &
                                   * MOMX(k+1,i,j)
             MOMX(k,i,j) = MOMX(k,i,j) - dt * ( eflx(k,i,j) - eflx(k-1,i,j) ) * rcdz_u(k,i,j)
          enddo
          MOMX(KE,i,j) = MOMX(KE,i,j) - dt * ( - eflx(KE-1,i,j) ) * rcdz_u(KE,i,j)

          !--- momentum y
          eflx(KS-1,i,j) = 0.25_RP * ( vterm(KS,i,j,iq ) + vterm(KS,i,j+1,iq ) ) &
                                   * ( QTRC (KS,i,j,iqa) + QTRC (KS,i,j+1,iqa) ) &
                                   * MOMY(KS,i,j)
          do k  = KS, KE-1
             eflx(k,i,j) = 0.25_RP * ( vterm(k+1,i,j,iq ) + vterm(k+1,i,j+1,iq ) ) &
                                   * ( QTRC (k+1,i,j,iqa) + QTRC (k+1,i,j+1,iqa) ) &
                                   * MOMY(k+1,i,j)
             MOMY(k,i,j) = MOMY(k,i,j) - dt * ( eflx(k,i,j) - eflx(k-1,i,j) ) * rcdz_v(k,i,j)
          enddo
          MOMY(KE,i,j) = MOMY(KE,i,j) - dt * ( - eflx(KE-1,i,j) ) * rcdz_v(KE,i,j)

       enddo
       enddo

    enddo ! falling (water mass & number) tracer

    !--- save previous value
    do iqa = 1, QA
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          rhoq(k,i,j,iqa) = QTRC(k,i,j,iqa) * DENS(k,i,j)
       enddo
       enddo
       enddo
    enddo ! all tracer

    do iq = 1, QA_MP-1
       iqa = QS_MP + iq

       !--- mass flux for each tracer, upwind with vel < 0
       do j = JS, JE
       do i = IS, IE
          do k  = KS-1, KE-1
             qflx(k,i,j,iq) = vterm(k+1,i,j,iq) * rhoq(k+1,i,j,iqa)
          enddo
          qflx(KE,i,j,iq) = 0.0_RP
       enddo
       enddo

       if ( TRACER_MASS(iqa) == 0.0_RP ) cycle

       !--- update total density
       do j  = JS, JE
       do i  = IS, IE
       do k  = KS, KE
          DENS(k,i,j) = DENS(k,i,j) - dt * ( qflx(k,i,j,iq) - qflx(k-1,i,j,iq) ) * rcdz(k,i,j)
       enddo
       enddo
       enddo

    enddo ! falling (water mass & number) tracer

    !--- update falling tracer
    do iq = 1, QA_MP-1
       iqa = QS_MP + iq
       do j  = JS, JE
       do i  = IS, IE
       do k  = KS, KE
          rhoq(k,i,j,iqa) = rhoq(k,i,j,iqa) - dt * ( qflx(k,i,j,iq) - qflx(k-1,i,j,iq) ) * rcdz(k,i,j)
       enddo
       enddo
       enddo
    enddo ! falling (water mass & number) tracer

    !--- update tracer ratio with updated total density)
    do iqa = 1, QA
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          QTRC(k,i,j,iqa) = rhoq(k,i,j,iqa) / DENS(k,i,j)
       enddo
       enddo
       enddo
    enddo ! all tracer

    call PROF_rapend  ('MP_Precipitation', 2)

    return
  end subroutine ATMOS_PHY_MP_precipitation_obsolute

end module scale_atmos_phy_mp_common
