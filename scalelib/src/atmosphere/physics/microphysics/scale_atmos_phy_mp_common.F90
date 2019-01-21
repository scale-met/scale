!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics - Common
!!
!! @par Description
!!          Common module for Cloud Microphysics
!!          Sedimentation/Precipitation and Saturation adjustment
!!
!! @author Team SCALE
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
  use scale_io
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
  public :: ATMOS_PHY_MP_precipitation_upwind
  public :: ATMOS_PHY_MP_precipitation_semilag
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
    !$omp        QV,QTRC,diffq_check)
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

    enddo
    enddo
    enddo


    diffq_min = minval( diffq_check(KS:KE,IS:IE,JS:JE) )

    if (       abs(limit_negative) > 0.0_RP         &
         .AND. abs(limit_negative) < abs(diffq_min) ) then
       LOG_ERROR("ATMOS_PHY_MP_negative_fixer",*) 'large negative is found. rank = ', PRC_myrank

       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          if (     abs(limit_negative) < abs(diffq_check(k,i,j)) &
              .OR. abs(QV(k,i,j)     ) < abs(diffq_check(k,i,j)) ) then
             LOG_ERROR_CONT(*) 'k,i,j,value(QHYD,QV) = ', k, i, j, diffq_check(k,i,j), QV(k,i,j)
          endif
       enddo
       enddo
       enddo
       LOG_ERROR_CONT(*) 'maximum negative hydrometeor ', diffq_min, ' < ', - abs(limit_negative)

       call PRC_abort
    endif

    !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
    !$omp private(i,j,k, &
    !$omp         diffq) &
    !$omp shared(KS,KE,IS,IE,JS,JE, &
    !$omp        QV,DENS)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
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

    call PROF_rapend('MP_filter', 3)

    return
  end subroutine ATMOS_PHY_MP_negative_fixer

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
!OCL SERIAL
  subroutine ATMOS_PHY_MP_precipitation_upwind( &
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

    real(RP) :: vtermh(KA)
    real(RP) :: qflx  (KA)
    real(RP) :: eflx  (KA)
    real(RP) :: RHOCP (KA)
    real(RP) :: RHOCV (KA)
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
       do k = KS, KE-1
          vtermh(k) = 0.5_RP * ( vterm(k+1,iq) + vterm(k,iq) )
       enddo
       vtermh(KS-1) = vterm(KS,iq)

       !--- mass flux for each tracer, upwind with vel < 0
       do k = KS-1, KE-1
          qflx(k) = vtermh(k) * RHOQ(k+1,iq)
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
  end subroutine ATMOS_PHY_MP_precipitation_upwind

  !-----------------------------------------------------------------------------
!OCL SERIAL
  subroutine ATMOS_PHY_MP_precipitation_semilag( &
       KA, KS, KE, QHA, QLA, QIA, &
       TEMP, vterm, FZ, FDZ, RCDZ, dt, &
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
    real(RP), intent(in) :: FZ   (KA)
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

    real(RP) :: vtermh(KA)
    real(RP) :: dvterm(KA)
    real(RP) :: cdz   (KA)
    real(RP) :: rfdz2 (KA)
    real(RP) :: dist  (KA)
    real(RP) :: Z_src
    real(RP) :: flx
    integer  :: k_src (KA)
    integer  :: k_dst

    integer  :: k, iq
    !---------------------------------------------------------------------------

    ! tracer/energy transport by falldown
    ! velocity is always negative

    mflx(:) = 0.0_RP
    sflx(:) = 0.0_RP
    qflx(:) = 0.0_RP
    eflx(:) = 0.0_RP

    do k = KS, KE
       RHOCP(k) = CPtot(k) * DENS(k)
       RHOCV(k) = CVtot(k) * DENS(k)
    end do

    do k = KS, KE
       CDZ(k) = 1.0_RP / RCDZ(k)
    end do
    do k = KS, KE-1
       rfdz2(k) = 1.0_RP / ( CDZ(k) + CDZ(k+1) )
    end do

    do iq = 1, QHA
       do k = KS, KE-1
          vtermh(k) = ( CDZ(k) * vterm(k+1,iq) + CDZ(k+1) * vterm(k,iq) ) * rfdz2(k)
       enddo
       vtermh(KS-1) = vterm(KS,iq)

       do k = KS, KE
          dvterm(k) = vtermh(k) - vtermh(k-1)
       enddo

       ! Movement distance of the cell wall by the fall
       ! the midpoint method (second-order Runge-Kutta)
       ! dz/dt = v(z + v dt/2) ~ v(z) + v dt/2 dv/dz + 1/2 (v dt/2)^2 d^2v/dz^2
       do k = KS, KE-1
          dist(k) = - vtermh(k)    * dt                                                                        &
                    + vtermh(k)    * dt**2 / 2.0_RP * ( dvterm(k+1)+dvterm(k) ) * rfdz2(k)                     &
                    - vtermh(k)**2 * dt**3 / 4.0_RP * ( dvterm(k+1)*RCDZ(k+1) - dvterm(k)*RCDZ(k) ) * rfdz2(k)
          dist(k) = max( dist(k), 0.0_RP )
       enddo
       dist(KS-1) = - vtermh(KS-1) * dt &
                    + vtermh(KS-1) * dt**2 / 2.0_RP * dvterm(KS)*RCDZ(KS)
       dist(KS-1) = max( dist(KS-1), 0.0_RP )

       ! wall cannot overtake
       do k = KE-2, KS-1, -1
          dist(k) = min( dist(k), dist(k+1) + CDZ(k+1) )
       end do

!        LOG_INFO_CONT(*) "distance", iq
!        do k = KA, 1, -1
!           LOG_INFO_CONT('(1x,I5,3F9.3,ES15.5)') k, dist(k), vtermh(k), vterm(k,iq), RHOQ(k,iq)
!        enddo

       ! search number of source cell
       do k_dst = KS-1, KE-1
          Z_src = FZ(k_dst) + dist(k_dst)

          k_src(k_dst) = k_dst
          do k = k_dst, KE-1
             if (       Z_src >  FZ(k  ) &
                  .AND. Z_src <= FZ(k+1) ) then
                k_src(k_dst) = k
             endif
          enddo
          if ( Z_src > FZ(KE) ) k_src(k_dst) = KE
       enddo

!        LOG_INFO_CONT(*) "seek", iq
!        do k = KA, 1, -1
!           LOG_INFO_CONT('(1x,2I5,2F9.3)') k, k_src(k), FZ(k), FZ(k)+dist(k)
!        enddo

       if ( iq > QLA ) then ! ice water
          CP = CP_ICE
          CV = CV_ICE
       else                 ! liquid water
          CP = CP_WATER
          CV = CV_WATER
       end if

       do k_dst = KS-1, KE-1
          do k = k_dst, k_src(k_dst)-1
             flx = RHOQ(k+1,iq) * CDZ(k+1) / dt               ! sum column mass rhoq*dz
             qflx(k_dst) = qflx(k_dst) - flx
             eflx(k_dst) = eflx(k_dst) - flx * TEMP(k+1) * CV ! internal energy flux
             dist(k_dst) = dist(k_dst) - CDZ(k+1)             ! residual
          enddo
          if ( k_src(k_dst) < KE ) then
             ! residual (simple upwind)
             flx = RHOQ(k_src(k_dst)+1,iq) * dist(k_dst) / dt
             qflx(k_dst) = qflx(k_dst) - flx                            ! sum column mass rhoq*dz
             eflx(k_dst) = eflx(k_dst) -flx * TEMP(k_src(k_dst)+1) * CV ! internal energy flux
          end if

          eflx(k_dst) = eflx(k_dst) + qflx(k) * FDZ(k_dst) * GRAV ! potential energy
       enddo

!        LOG_INFO_CONT(*) "flux", iq
!        do k = KA, 1, -1
!           LOG_INFO_CONT('(1x,2I5,F9.3,2ES15.5)') k, k_src(k), dist(k), qflx(k), vtermh(k)*RHOQ(k+1,iq)
!        enddo

       !--- update falling tracer
       do k  = KS, KE
          rhoq(k,iq) = rhoq(k,iq) - dt * ( qflx(k) - qflx(k-1) ) * RCDZ(k)
       enddo ! falling (water mass & number) tracer

!        LOG_INFO_CONT(*) "tendency", iq
!        do k = KA, 1, -1
!           LOG_INFO_CONT('(1x,I5,ES15.5)') k, - dt * ( qflx(k) - qflx(k-1) ) * RCDZ(k)
!        enddo

       ! QTRC(iq; iq>QLA+QLI) is not mass tracer, such as number density
       if ( iq > QLA + QIA ) cycle

       do k = KS-1, KE-1
          mflx(k) = mflx(k) + qflx(k)
       end do

       if ( iq > QLA ) then ! ice water
          sflx(2) = sflx(2) + qflx(KS-1)
       else                 ! liquid water
          sflx(1) = sflx(1) + qflx(KS-1)
       end if

       !--- update density
       do k = KS, KE
          dDENS = - ( qflx(k) - qflx(k-1) ) * RCDZ(k) * dt
          RHOCP(k) = RHOCP(k) + CP * dDENS
          RHOCV(k) = RHOCV(k) + CV * dDENS
          DENS(k) = DENS(k) + dDENS
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
  end subroutine ATMOS_PHY_MP_precipitation_semilag

  !-----------------------------------------------------------------------------
!OCL SERIAL
  subroutine ATMOS_PHY_MP_precipitation_momentum( &
       KA, KS, KE, &
       DENS, MOMZ, U, V, mflx, &
       RCDZ, RFDZ,             &
       MOMZ_t, RHOU_t, RHOV_t  )
    implicit none

    integer,  intent(in)  :: KA, KS, KE
    real(RP), intent(in)  :: DENS  (KA)
    real(RP), intent(in)  :: MOMZ  (KA)
    real(RP), intent(in)  :: U     (KA)
    real(RP), intent(in)  :: V     (KA)
    real(RP), intent(in)  :: mflx  (KA)
    real(RP), intent(in)  :: RCDZ  (KA)
    real(RP), intent(in)  :: RFDZ  (KA)
    real(RP), intent(out) :: MOMZ_t(KA)
    real(RP), intent(out) :: RHOU_t(KA)
    real(RP), intent(out) :: RHOV_t(KA)

    real(RP) :: flx(KA)

    integer  :: k
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

end module scale_atmos_phy_mp_common
