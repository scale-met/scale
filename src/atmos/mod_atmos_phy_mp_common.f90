!-------------------------------------------------------------------------------
!> module Atmosphere / Physics Cloud Microphysics
!!
!! @par Description
!!          Common module for Cloud Microphysics
!!          Sedimentation/Precipitation and Saturation adjustment
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2012-12-23 (H.Yashiro)  [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_phy_mp_common
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_stdio, only: &
     IO_FID_LOG,  &
     IO_L
  use mod_time, only: &
     TIME_rapstart, &
     TIME_rapend
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

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include "inc_precision.h"
  include 'inc_index.h'
  include 'inc_tracer.h'

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  public :: moist_conversion_liq
  public :: moist_conversion_all

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Negative fixer
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_negative_fixer( &
       DENS, &
       RHOT, &
       QTRC  )
    implicit none

    real(RP), intent(inout) :: DENS(KA,IA,JA)
    real(RP), intent(inout) :: RHOT(KA,IA,JA)
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA)

    real(RP) :: diffq(KA)

    integer :: i, j, iq
    !---------------------------------------------------------------------------

    do j = 1, JA
    do i = 1, IA
       diffq(:) = 0.D0

       ! total hydrometeor (before correction)
       do iq = QQS+1, QQE
          diffq(:) = diffq(:) + QTRC(:,i,j,iq)
       enddo

       ! remove negative value of hydrometeors (mass)
       do iq = QQS+1, QQE
          QTRC(:,i,j,iq) = max( QTRC(:,i,j,iq), 0.0_RP )
       enddo

       ! difference between before and after correction
       do iq = QQS+1, QQE
          diffq(:) = diffq(:) - QTRC(:,i,j,iq)
       enddo

       ! Compensate for the lack of hydrometeors by the water vapor
       QTRC(:,i,j,I_QV) = QTRC(:,i,j,I_QV) + diffq(:)

       ! TODO: We have to consider energy conservation (but very small)

       ! remove negative value of water vapor (mass)
       diffq(:) = QTRC(:,i,j,I_QV)
       QTRC(:,i,j,I_QV) = max( QTRC(:,i,j,I_QV), 0.0_RP )
       diffq(:) = diffq(:) - QTRC(:,i,j,I_QV)

       ! Apply correction to total density
       ! TODO: We have to consider energy conservation (but very small)
       DENS(:,i,j) = DENS(:,i,j) * ( 1.0_RP + diffq(:) )
       RHOT(:,i,j) = RHOT(:,i,j) * ( 1.0_RP + diffq(:) )

    enddo
    enddo

    return
  end subroutine ATMOS_PHY_MP_negative_fixer

  !-----------------------------------------------------------------------------
  !> Saturation adjustment
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_saturation_adjustment( &
       RHOT0, &
       QTRC0, &
       DENS0  )
    use mod_const, only: &
       LHV00  => CONST_LH00, &
       LHF00  => CONST_LHF00
    use mod_atmos_thermodyn, only: &
       THERMODYN_qd        => ATMOS_THERMODYN_qd,       &
       THERMODYN_cv        => ATMOS_THERMODYN_cv,       &
       THERMODYN_rhot      => ATMOS_THERMODYN_rhot,     &
       THERMODYN_temp_pres => ATMOS_THERMODYN_temp_pres
    use mod_atmos_saturation, only: &
       SATURATION_dens2qsat_liq => ATMOS_SATURATION_dens2qsat_liq, &
       SATURATION_dens2qsat_all => ATMOS_SATURATION_dens2qsat_all
    implicit none

    real(RP), intent(inout) :: RHOT0(KA,IA,JA)
    real(RP), intent(inout) :: QTRC0(KA,IA,JA,QA)
    real(RP), intent(in)    :: DENS0(KA,IA,JA)

    ! working
    real(RP) :: dens      ! density [kg/m3]
    real(RP) :: rhot      ! density * potential temperature [K*kg/m3]
    real(RP) :: q(QA)     ! tracer Q [kg/kg]
    real(RP) :: temp      ! temperature [K]
    real(RP) :: pres      ! pressure [Pa]
    real(RP) :: rhoe      ! internal energy

    real(RP) :: qsum      ! QV+QC+QI
    real(RP) :: qdry
    real(RP) :: ein_moist ! moist internal energy

    real(RP) :: qsat      ! saturated water vapor
    real(RP) :: CVtot

    integer :: conversion_count

    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    call TIME_rapstart('MP_saturation_adjustment')

    conversion_count = 0

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE

       ! store to 1D work
       dens = DENS0(k,i,j)
       rhot = RHOT0(k,i,j)
       do iq = QQS, QQE
          q(iq) = QTRC0(k,i,j,iq)
       enddo

       call THERMODYN_temp_pres( temp, & ! [OUT]
                                 pres, & ! [OUT]
                                 dens, & ! [IN]
                                 rhot, & ! [IN]
                                 q(:)  ) ! [IN]


       !##### Start Main #####

       ! Turn QC & QI into QV with consistency of moist internal energy
       call THERMODYN_qd( qdry,  q(:) ) ! qdry dont change through the process
       call THERMODYN_cv( CVtot, q(:), qdry )

       if ( I_QI <= 0 ) then ! warm rain

          ein_moist = temp * CVtot + q(I_QV) * LHV00
          qsum    = q(I_QV) + q(I_QC)
          q(I_QV) = qsum
          q(I_QC) = 0.0_RP

          call THERMODYN_cv( CVtot, q(:), qdry )

          ! new temperature (after QC evaporation)
          temp = ( ein_moist - qsum * LHV00 ) / CVtot

          call SATURATION_dens2qsat_liq( qsat, temp, dens )

          ! saturation adjustment
          if ( qsum > qsat ) then
             conversion_count = conversion_count + 1

             call moist_conversion_liq( temp,     & ! [INOUT]
                                        q(:),     & ! [INOUT]
                                        dens,     & ! [IN]
                                        qsum,     & ! [IN]
                                        qdry,     & ! [IN]
                                        ein_moist ) ! [IN]
          endif

       else ! cold rain

          ein_moist = temp * CVtot + q(I_QV) * LHV00 - q(I_QI) * LHF00
          qsum    = q(I_QV) + q(I_QC) + q(I_QI)
          q(I_QV) = qsum
          q(I_QC) = 0.0_RP
          q(I_QI) = 0.0_RP

          call THERMODYN_cv( CVtot, q(:), qdry )

          ! new temperature (after QC & QI evaporation)
          temp = ( ein_moist - qsum * LHV00 ) / CVtot

          call SATURATION_dens2qsat_all( qsat, temp, dens )

          ! saturation adjustment
          if ( qsum > qsat ) then
             conversion_count = conversion_count + 1

             call moist_conversion_all( temp,     & ! [INOUT]
                                        q(:),     & ! [INOUT]
                                        dens,     & ! [IN]
                                        qsum,     & ! [IN]
                                        qdry,     & ! [IN]
                                        ein_moist ) ! [IN]
          endif

       endif

       call THERMODYN_cv( CVtot, q(:), qdry )

       rhoe = dens * temp * CVtot

       !##### End Main #####


       call THERMODYN_rhot( rhot, & ! [OUT]
                            rhoe, & ! [IN]
                            q(:)  ) ! [IN]

       RHOT0(k,i,j) = rhot
       do iq = QQS, QQE
          QTRC0(k,i,j,iq) = q(iq)
       enddo

    enddo
    enddo
    enddo

    if( IO_L ) write(IO_FID_LOG,*) '*** # of conversion point = ', conversion_count

    call TIME_rapend  ('MP_saturation_adjustment')

    return
  end subroutine ATMOS_PHY_MP_saturation_adjustment

  !-----------------------------------------------------------------------------
  !> Iterative moist conversion for warm rain
  !-----------------------------------------------------------------------------
  subroutine moist_conversion_liq( &
       temp,      &
       q,         &
       dens,      &
       qsum,      &
       qdry,      &
       ein_moist0 )
    use mod_const, only : &
       LHV00  => CONST_LH00
    use mod_process, only: &
       PRC_MPIstop
    use mod_atmos_thermodyn, only: &
       THERMODYN_cv => ATMOS_THERMODYN_cv, &
       CVw => AQ_CV
    use mod_atmos_saturation, only: &
       SATURATION_dens2qsat_liq => ATMOS_SATURATION_dens2qsat_liq, &
       SATURATION_alpha         => ATMOS_SATURATION_alpha,         &
       SATURATION_dalphadT      => ATMOS_SATURATION_dalphadT,      &
       CVovR_liq, &
       LovR_liq
    implicit none

    real(RP), intent(inout) :: temp
    real(RP), intent(inout) :: q(QA)
    real(RP), intent(in)    :: dens
    real(RP), intent(in)    :: qsum
    real(RP), intent(in)    :: qdry
    real(RP), intent(in)    :: ein_moist0

    ! working
    real(RP) :: dtemp

    real(RP) :: qsatl
    real(RP) :: CVtot
    real(RP) :: ein_moist ! moist internal energy

    ! d(X)/dT
    real(RP) :: dqsatl_dT
    real(RP) :: dqc_dT
    real(RP) :: dCVtot_dT
    real(RP) :: dein_moist_dT

    integer  :: itelim = 10
    real(RP) :: dtemp_criteria = 1.E-6_RP

    integer :: ite
    !---------------------------------------------------------------------------

    do ite = 1, itelim

       call SATURATION_dens2qsat_liq( qsatl, temp, dens )

       ! Separation
       q(I_QV) = qsatl
       q(I_QC) = qsum-qsatl

       call THERMODYN_cv( CVtot, q(:), qdry )

       ein_moist = temp * CVtot + qsatl * LHV00

       ! dX/dT
       dqsatl_dT = ( LovR_liq / ( temp*temp ) + CVovR_liq / temp ) * qsatl

       dqc_dT = - dqsatl_dT

       dCVtot_dT = dqsatl_dT * CVw(I_QV) &
                 + dqc_dT    * CVw(I_QC)

       dein_moist_dT = temp * dCVtot_dT + CVtot + dqsatl_dT * LHV00

       dtemp = ( ein_moist - ein_moist0 ) / dein_moist_dT
       temp  = temp - dtemp

       if ( abs(dtemp) < dtemp_criteria ) exit

    enddo

    if ( ite > itelim ) then
       if( IO_L ) write(IO_FID_LOG,*) 'xxx [moist_conversion] not converged! dtemp=', dtemp
       call PRC_MPIstop
    endif

    return
  end subroutine moist_conversion_liq

  !-----------------------------------------------------------------------------
  !> Iterative moist conversion (liquid/ice mixture)
  !-----------------------------------------------------------------------------
  subroutine moist_conversion_all( &
       temp,      &
       q,         &
       dens,      &
       qsum,      &
       qdry,      &
       ein_moist0 )
    use mod_const, only : &
       LHV00  => CONST_LH00,  &
       LHF00  => CONST_LHF00
    use mod_process, only: &
       PRC_MPIstop
    use mod_atmos_thermodyn, only: &
       THERMODYN_cv => ATMOS_THERMODYN_cv, &
       CVw => AQ_CV
    use mod_atmos_saturation, only: &
       SATURATION_dens2qsat_all => ATMOS_SATURATION_dens2qsat_all, &
       SATURATION_dens2qsat_liq => ATMOS_SATURATION_dens2qsat_liq, &
       SATURATION_dens2qsat_ice => ATMOS_SATURATION_dens2qsat_ice, &
       SATURATION_alpha         => ATMOS_SATURATION_alpha,         &
       SATURATION_dalphadT      => ATMOS_SATURATION_dalphadT,      &
       CVovR_liq, &
       CVovR_ice, &
       LovR_liq,  &
       LovR_ice
    implicit none

    real(RP), intent(inout) :: temp
    real(RP), intent(inout) :: q(QA)
    real(RP), intent(in)    :: dens
    real(RP), intent(in)    :: qsum
    real(RP), intent(in)    :: qdry
    real(RP), intent(in)    :: ein_moist0

    real(RP) :: alpha
    real(RP) :: qsat, qsatl, qsati
    real(RP) :: CVtot
    real(RP) :: ein_moist ! moist internal energy

    ! d(X)/d(T)
    real(RP) :: dalpha_dT
    real(RP) :: dqsat_dT, dqsatl_dT, dqsati_dT
    real(RP) :: dCVtot_dT
    real(RP) :: dein_moist_dT

    real(RP) :: dqc_dT, dqi_dT
    real(RP) :: dtemp

    integer  :: itelim = 10
    real(RP) :: dtemp_criteria = 1.E-6_RP

    integer :: ite
    !---------------------------------------------------------------------------

    do ite = 1, itelim

       ! liquid/ice separation factor
       call SATURATION_alpha( alpha, temp )
       ! Saturation
       call SATURATION_dens2qsat_all( qsat,  temp, dens )
       call SATURATION_dens2qsat_liq( qsatl, temp, dens )
       call SATURATION_dens2qsat_ice( qsati, temp, dens )

       ! Separation
       q(I_QV) = qsat
       q(I_QC) = ( qsum-qsat ) * (        alpha )
       q(I_QI) = ( qsum-qsat ) * ( 1.0_RP-alpha )

       call THERMODYN_cv( CVtot, q(:), qdry )

       ein_moist = temp * CVtot + qsat * LHV00 - q(I_QI) * LHF00

       ! dX/dT
       call SATURATION_dalphadT( dalpha_dT, temp )

       dqsatl_dT = ( LovR_liq / ( temp*temp ) + CVovR_liq / temp ) * qsatl
       dqsati_dT = ( LovR_ice / ( temp*temp ) + CVovR_ice / temp ) * qsati

       dqsat_dT  = qsatl * dalpha_dT + dqsatl_dT * (        alpha ) &
                 - qsati * dalpha_dT + dqsati_dT * ( 1.0_RP-alpha )

       dqc_dT =  ( qsum-qsat ) * dalpha_dT - dqsat_dT * (        alpha )
       dqi_dT = -( qsum-qsat ) * dalpha_dT - dqsat_dT * ( 1.0_RP-alpha )

       dCVtot_dT = dqsat_dT * CVw(I_QV) &
                 + dqc_dT   * CVw(I_QC) &
                 + dqi_dT   * CVw(I_QI)

       dein_moist_dT = temp * dCVtot_dT + CVtot + dqsat_dT * LHV00 - dqi_dT * LHF00

       dtemp = ( ein_moist - ein_moist0 ) / dein_moist_dT
       temp  = temp - dtemp

       if ( abs(dtemp) < dtemp_criteria ) exit

    enddo

    if ( ite > itelim ) then
       if( IO_L ) write(IO_FID_LOG,*) 'xxx [moist_conversion] not converged! dtemp=', dtemp
       call PRC_MPIstop
    endif

    return
  end subroutine moist_conversion_all

  !-----------------------------------------------------------------------------
  !> precipitation transport
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_precipitation( &
       flux_rain, &
       flux_snow, &
       DENS,      &
       MOMZ,      &
       MOMX,      &
       MOMY,      &
       RHOE,      &
       QTRC,      &
       vterm,     &
       temp       )
    use mod_const, only : &
       GRAV  => CONST_GRAV
    use mod_time, only: &
       dt => TIME_DTSEC_ATMOS_PHY_MP
    use mod_grid, only: &
       CZ   => GRID_CZ,   &
       FDZ  => GRID_FDZ,  &
       RCDZ => GRID_RCDZ, &
       RFDZ => GRID_RFDZ
    use mod_atmos_thermodyn, only: &
       CVw => AQ_CV
    implicit none

    real(RP), intent(out)   :: flux_rain(KA,IA,JA)
    real(RP), intent(out)   :: flux_snow(KA,IA,JA)
    real(RP), intent(inout) :: DENS     (KA,IA,JA)
    real(RP), intent(inout) :: MOMZ     (KA,IA,JA)
    real(RP), intent(inout) :: MOMX     (KA,IA,JA)
    real(RP), intent(inout) :: MOMY     (KA,IA,JA)
    real(RP), intent(inout) :: RHOE     (KA,IA,JA)
    real(RP), intent(inout) :: QTRC     (KA,IA,JA,QA)
    real(RP), intent(in)    :: vterm    (KA,IA,JA,QA) ! terminal velocity of cloud mass
    real(RP), intent(in)    :: temp     (KA,IA,JA)

    real(RP) :: rhoq(KA,QA) ! rho * q before precipitation
    real(RP) :: qflx(KA,QA)
    real(RP) :: eflx(KA)

    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    call TIME_rapstart('MP_precipitation')

    flux_rain(:,:,:) = 0.0_RP
    flux_snow(:,:,:) = 0.0_RP

    ! tracer/energy transport by falldown
    ! 1st order upwind, forward euler, velocity is always negative
    do j  = JS, JE
    do i  = IS, IE

       eflx(KE) = 0.0_RP

       !--- mass flux for each mass tracer, upwind with vel < 0
       do iq = I_QC, QA
          do k  = KS-1, KE-1
             rhoq(k,iq) = DENS(k,i,j) * QTRC(k,i,j,iq)
             qflx(k,iq) = vterm(k+1,i,j,iq) * DENS(k+1,i,j) * QTRC(k+1,i,j,iq)
          enddo
          rhoq(KE,iq) = DENS(KE,i,j) * QTRC(KE,i,j,iq)
          qflx(KE,iq) = 0.0_RP
       enddo

       !--- lowermost flux is saved for land process
       do k  = KS-1, KE
          if ( QWS > 0 ) then
             do iq = QWS, QWE
                flux_rain(k,i,j) = flux_rain(k,i,j) - qflx(k,iq)
             enddo
          endif
          if ( QIS > 0 ) then
             do iq = QIS, QIE
                flux_snow(k,i,j) = flux_snow(k,i,j) - qflx(k,iq)
             enddo
          endif
       enddo

       do iq = I_QC, QQE ! mass tracer only

          !--- internal energy
          do k  = KS-1, KE-1
             eflx(k) = qflx(k,iq) * temp(k+1,i,j) * CVw(iq)
          enddo
          do k  = KS, KE
             RHOE(k,i,j) = RHOE(k,i,j) - dt * ( eflx(k) - eflx(k-1) ) * RCDZ(k)
          enddo

          !--- potential energy
          do k  = KS, KE-1
             eflx(k) = qflx(k,iq) * GRAV * FDZ(k)
          enddo
          eflx(KS-1) = qflx(KS-1,iq) * GRAV * CZ(KS)
          do k  = KS, KE
             RHOE(k,i,j) = RHOE(k,i,j) - dt * ( eflx(k) - eflx(k-1) ) * RCDZ(k)
          enddo

          !--- momentum z (half level)
          do k  = KS-1, KE-2
             eflx(k) = 0.25_RP * ( vterm(k+1,i,j,iq) + vterm(k,i,j,iq) ) &
                               * ( QTRC (k+1,i,j,iq) + QTRC (k,i,j,iq) ) &
                               * MOMZ(k,i,j)
          enddo
          do k  = KS, KE-1
             MOMZ(k,i,j) = MOMZ(k,i,j) - dt * ( eflx(k+1) - eflx(k) ) * RFDZ(k)
          enddo

          !--- momentum x
          do k  = KS-1, KE-1
             eflx(k) = 0.25_RP * ( vterm(k+1,i,j,iq) + vterm(k+1,i+1,j,iq) ) &
                               * ( QTRC (k+1,i,j,iq) + QTRC (k+1,i+1,j,iq) ) &
                               * MOMX(k+1,i,j)
          enddo
          do k  = KS, KE
             MOMX(k,i,j) = MOMX(k,i,j) - dt * ( eflx(k) - eflx(k-1) ) * RCDZ(k)
          enddo

          !--- momentum y
          do k  = KS-1, KE-1
             eflx(k) = 0.25_RP * ( vterm(k+1,i,j,iq) + vterm(k+1,i,j+1,iq) ) &
                               * ( QTRC (k+1,i,j,iq) + QTRC (k+1,i,j+1,iq) ) &
                               * MOMY(k+1,i,j)
          enddo
          do k  = KS, KE
             MOMY(k,i,j) = MOMY(k,i,j) - dt * ( eflx(k) - eflx(k-1) ) * RCDZ(k)
          enddo

          !--- update total density
          do k  = KS, KE
             DENS(k,i,j) = DENS(k,i,j) - dt * ( qflx(k,iq) - qflx(k-1,iq) ) * RCDZ(k)
          enddo

       enddo ! mass tracer loop

       !--- update tracer
       do iq = I_QC, QA
       do k  = KS, KE
          QTRC(k,i,j,iq) = ( rhoq(k,iq) - dt * ( qflx(k,iq) - qflx(k-1,iq) ) * RCDZ(k) ) &
                         / DENS(k,i,j)
       enddo 
       enddo

    enddo ! I loop
    enddo ! J loop

    call TIME_rapend  ('MP_precipitation')

    return
  end subroutine ATMOS_PHY_MP_precipitation

end module mod_atmos_phy_mp_common 
