!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics
!!
!! @par Description
!!          Cloud Microphysics by Kessler-type parametarization
!!          Reference: Kessler(1969)
!!                     Klemp and Wilhelmson(1978)
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2012-01-14 (Y.Miyamoto) [new]
!! @li      2012-03-23 (H.Yashiro)  [mod] Explicit index parameter inclusion
!! @li      2012-12-23 (H.Yashiro)  [mod] Reconstruction
!!
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module scale_atmos_phy_mp_kessler
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index

  use scale_tracer_kessler
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_MP_kessler_setup
  public :: ATMOS_PHY_MP_kessler
  public :: ATMOS_PHY_MP_kessler_CloudFraction
  public :: ATMOS_PHY_MP_kessler_EffectiveRadius
  public :: ATMOS_PHY_MP_kessler_Mixingratio

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, allocatable :: vterm(:,:,:,:) ! terminal velocity of each tracer [m/s]

  real(RP), public, target :: ATMOS_PHY_MP_DENS(MP_QA) ! hydrometeor density [kg/m3]=[g/L]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: MP_kessler
  private :: MP_kessler_vterm

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical :: MP_doreport_tendency = .false. ! report tendency of each process?
  logical :: MP_donegative_fixer  = .true. ! apply negative fixer?

  real(RP), allocatable :: factor_vterm(:) ! collection factor for terminal velocity of QR

  logical :: first = .true.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup Cloud Microphysics
  subroutine ATMOS_PHY_MP_kessler_setup( MP_TYPE )
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       CONST_DWATR
    implicit none
    character(len=H_SHORT), intent(in) :: MP_TYPE

    NAMELIST / PARAM_ATMOS_PHY_MP / &
       MP_doreport_tendency, &
       MP_donegative_fixer

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Cloud Microphisics]/Categ[ATMOS]'
    if( IO_L ) write(IO_FID_LOG,*) '*** KESSLER-type 1-moment bulk 3 category'

    allocate( vterm(KA,IA,JA,QA) )
    allocate( factor_vterm(KA) )


    if ( MP_TYPE /= 'KESSLER' ) then
       if ( IO_L ) write(IO_FID_LOG,*) 'xxx ATMOS_PHY_MP_TYPE is not KESSLER. Check!'
       call PRC_MPIstop
    endif

    if (      I_QV <= 0 &
         .OR. I_QC <= 0 &
         .OR. I_QR <= 0 ) then
       if ( IO_L ) write(IO_FID_LOG,*) 'xxx KESSLER needs QV, QC, QR tracer. Check!'
       call PRC_MPIstop
    endif

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_MP,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_MP. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_MP)

    vterm(:,:,:,:) = 0.0_RP

    ATMOS_PHY_MP_DENS(I_mp_QC) = CONST_DWATR
    ATMOS_PHY_MP_DENS(I_mp_QR) = CONST_DWATR

    return
  end subroutine ATMOS_PHY_MP_kessler_setup

  !-----------------------------------------------------------------------------
  !> Cloud Microphysics
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_kessler( &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       RHOT, &
       QTRC  )
    use scale_time, only: &
       dt => TIME_DTSEC_ATMOS_PHY_MP
    use scale_history, only: &
       HIST_in
    use scale_comm, only: &
       COMM_horizontal_mean
    use scale_atmos_phy_mp_common, only: &
       MP_negative_fixer        => ATMOS_PHY_MP_negative_fixer,       &
       MP_precipitation         => ATMOS_PHY_MP_precipitation,        &
       MP_saturation_adjustment => ATMOS_PHY_MP_saturation_adjustment
    use scale_atmos_thermodyn, only: &
       THERMODYN_rhoe        => ATMOS_THERMODYN_rhoe,       &
       THERMODYN_rhot        => ATMOS_THERMODYN_rhot,       &
       THERMODYN_temp_pres_E => ATMOS_THERMODYN_temp_pres_E
    use scale_tracer, only: &
       QAD => QA
    implicit none

    real(RP), intent(inout) :: DENS(KA,IA,JA)
    real(RP), intent(inout) :: MOMZ(KA,IA,JA)
    real(RP), intent(inout) :: MOMX(KA,IA,JA)
    real(RP), intent(inout) :: MOMY(KA,IA,JA)
    real(RP), intent(inout) :: RHOT(KA,IA,JA)
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QAD)

    real(RP) :: RHOE_t(KA,IA,JA)
    real(RP) :: QTRC_t(KA,IA,JA,QA)
    real(RP) :: RHOE  (KA,IA,JA)
    real(RP) :: temp  (KA,IA,JA)
    real(RP) :: pres  (KA,IA,JA)

    real(RP) :: flux_tot (KA,IA,JA)
    real(RP) :: flux_rain(KA,IA,JA)
    real(RP) :: flux_snow(KA,IA,JA)

    real(RP) :: rho_prof(KA) ! averaged profile of rho

    integer :: k
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Microphysics(kessler)'

    if ( first ) then
       ! Calculate collection factor for terminal velocity of QR
       call COMM_horizontal_mean( rho_prof(:), DENS(:,:,:) )
       rho_prof(:) = rho_prof(:) * 1.E-3_RP ! [kg/m3]->[g/cc]

       do k = KS, KE
          factor_vterm(k) = sqrt( rho_prof(KS)/rho_prof(k) )
       enddo
       first = .false.
    end if



    if ( MP_donegative_fixer ) then
       call MP_negative_fixer( DENS(:,:,:),  & ! [INOUT]
                               RHOT(:,:,:),  & ! [INOUT]
                               QTRC(:,:,:,:) ) ! [INOUT]
    endif

    call THERMODYN_rhoe( RHOE(:,:,:),  & ! [OUT]
                         RHOT(:,:,:),  & ! [IN]
                         QTRC(:,:,:,:) ) ! [IN]

    !##### MP Main #####
    RHOE_t(:,:,:)   = 0.0_RP
    QTRC_t(:,:,:,:) = 0.0_RP

    call MP_kessler( RHOE_t(:,:,:),   & ! [INOUT]
                     QTRC_t(:,:,:,:), & ! [INOUT]
                     RHOE  (:,:,:),   & ! [INOUT]
                     QTRC  (:,:,:,:), & ! [INOUT]
                     DENS  (:,:,:)    ) ! [IN]

    call MP_kessler_vterm( vterm(:,:,:,:), & ! [OUT]
                           DENS (:,:,:),   & ! [IN]
                           QTRC (:,:,:,:)  ) ! [IN]

    call THERMODYN_temp_pres_E( temp(:,:,:),  & ! [OUT]
                                pres(:,:,:),  & ! [OUT]
                                DENS(:,:,:),  & ! [IN]
                                RHOE(:,:,:),  & ! [IN]
                                QTRC(:,:,:,:) ) ! [IN]

    call MP_precipitation( flux_rain(:,:,:), & ! [OUT]
                           flux_snow(:,:,:), & ! [OUT]
                           DENS (:,:,:),     & ! [INOUT]
                           MOMZ (:,:,:),     & ! [INOUT]
                           MOMX (:,:,:),     & ! [INOUT]
                           MOMY (:,:,:),     & ! [INOUT]
                           RHOE (:,:,:),     & ! [INOUT]
                           QTRC (:,:,:,:),   & ! [INOUT]
                           vterm(:,:,:,:),   & ! [IN]
                           temp (:,:,:),     & ! [IN]
                           dt                )

    call MP_saturation_adjustment( RHOE_t(:,:,:),   & ! [INOUT]
                                   QTRC_t(:,:,:,:), & ! [INOUT]
                                   RHOE  (:,:,:),   & ! [INOUT]
                                   QTRC  (:,:,:,:), & ! [INOUT]
                                   DENS  (:,:,:)    ) ! [IN]

    if ( MP_doreport_tendency ) then
       call HIST_in( QTRC_t(:,:,:,I_QV), 'QV_t_mp', 'tendency of QV in mp', 'kg/kg/s', dt )
       call HIST_in( QTRC_t(:,:,:,I_QC), 'QC_t_mp', 'tendency of QC in mp', 'kg/kg/s', dt )
       call HIST_in( QTRC_t(:,:,:,I_QR), 'QR_t_mp', 'tendency of QR in mp', 'kg/kg/s', dt )

       call HIST_in( RHOE_t(:,:,:), 'RHOE_t_mp', 'tendency of rhoe in mp', 'J/m3/s', dt )

       call HIST_in( vterm(:,:,:,I_QR), 'Vterm_QR', 'terminal velocity of QR', 'm/s', dt )
    endif

    !##### END MP Main #####

    call THERMODYN_rhot( RHOT(:,:,:),  & ! [OUT]
                         RHOE(:,:,:),  & ! [IN]
                         QTRC(:,:,:,:) ) ! [IN]

    if ( MP_donegative_fixer ) then
       call MP_negative_fixer( DENS(:,:,:),  & ! [INOUT]
                               RHOT(:,:,:),  & ! [INOUT]
                               QTRC(:,:,:,:) ) ! [INOUT]
    endif

    flux_tot(:,:,:) = flux_rain(:,:,:) + flux_snow(:,:,:)
    call HIST_in( flux_rain(KS-1,:,:), 'RAIN', 'surface rain rate', 'kg/m2/s', dt)
    call HIST_in( flux_snow(KS-1,:,:), 'SNOW', 'surface snow rate', 'kg/m2/s', dt)
    call HIST_in( flux_tot (KS-1,:,:), 'PREC', 'surface precipitation rate', 'kg/m2/s', dt)

    return
  end subroutine ATMOS_PHY_MP_kessler

  !-----------------------------------------------------------------------------
  !> Kessler-type warm rain microphysics
  !-----------------------------------------------------------------------------
  subroutine MP_kessler( &
       RHOE_t, &
       QTRC_t, &
       RHOE0,  &
       QTRC0,  &
       DENS0   )
    use scale_const, only: &
       LHV00 => CONST_LH00
    use scale_time, only: &
       dt => TIME_DTSEC_ATMOS_PHY_MP
    use scale_atmos_thermodyn, only: &
       THERMODYN_temp_pres_E => ATMOS_THERMODYN_temp_pres_E
    use scale_atmos_saturation, only: &
       SATURATION_dens2qsat_liq => ATMOS_SATURATION_dens2qsat_liq
    implicit none

    real(RP), intent(inout) :: RHOE_t(KA,IA,JA)    ! tendency rhoe             [J/m3/s]
    real(RP), intent(inout) :: QTRC_t(KA,IA,JA,QA) ! tendency tracer           [kg/kg/s]
    real(RP), intent(inout) :: RHOE0 (KA,IA,JA)    ! density * internal energy [J/m3]
    real(RP), intent(inout) :: QTRC0 (KA,IA,JA,QA) ! mass concentration        [kg/kg]
    real(RP), intent(in)    :: DENS0 (KA,IA,JA)    ! density                   [kg/m3]

    ! working
    real(RP) :: QSAT (KA,IA,JA) ! saturated water vapor [kg/kg]
    real(RP) :: TEMP0(KA,IA,JA) ! temperature           [K]
    real(RP) :: PRES0(KA,IA,JA) ! pressure              [Pa]

    real(RP) :: dens (KA)
    real(RP) :: rhoe (KA)
    real(RP) :: temp (KA)
    real(RP) :: pres (KA)
    real(RP) :: qv   (KA)
    real(RP) :: qc   (KA)
    real(RP) :: qr   (KA)
    real(RP) :: qsatl(KA)
    real(RP) :: Sliq (KA)

    ! tendency
    real(RP) :: dq_evap(KA) ! tendency q (evaporation)
    real(RP) :: dq_auto(KA) ! tendency q (autoconversion)
    real(RP) :: dq_accr(KA) ! tendency q (accretion)

    real(RP) :: dqv, dqc, dqr
    real(RP) :: vent_factor
    real(RP) :: rdt

    integer :: k, i, j
    !---------------------------------------------------------------------------

    call PROF_rapstart('MP_kessler')

    rdt = 1.0_RP / dt

    call THERMODYN_temp_pres_E( TEMP0(:,:,:),  & ! [OUT]
                                PRES0(:,:,:),  & ! [OUT]
                                DENS0(:,:,:),  & ! [IN]
                                RHOE0(:,:,:),  & ! [IN]
                                QTRC0(:,:,:,:) ) ! [IN]

    call SATURATION_dens2qsat_liq( QSAT (:,:,:), & ! [OUT]
                                   TEMP0(:,:,:), & ! [IN]
                                   DENS0(:,:,:)  ) ! [IN]

    !$omp parallel do private(i,j,k,dens,rhoe,temp,pres,qv,qc,qr,qsatl,Sliq,dq_evap,dq_auto,dq_accr,vent_factor,dqc,dqr,dqv) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE

       ! store to work
       dens (KS:KE) = DENS0(KS:KE,i,j)
       rhoe (KS:KE) = RHOE0(KS:KE,i,j)
       temp (KS:KE) = TEMP0(KS:KE,i,j)
       pres (KS:KE) = PRES0(KS:KE,i,j)
       qv   (KS:KE) = QTRC0(KS:KE,i,j,I_QV)
       qc   (KS:KE) = QTRC0(KS:KE,i,j,I_QC)
       qr   (KS:KE) = QTRC0(KS:KE,i,j,I_QR)
       qsatl(KS:KE) = QSAT (KS:KE,i,j)
       Sliq (KS:KE) = qv(KS:KE) / qsatl(KS:KE)

       ! Auto-conversion (QC->QR)
       do k = KS, KE
          dq_auto(k) = 1.E-3_RP * max( qc(k)-1.E-3_RP, 0.0_RP )
       enddo

       ! Accretion (QC->QR)
       do k = KS, KE
          dq_accr(k) = 2.2_RP * qc(k) * qr(k)**0.875_RP
       enddo

       ! Evaporation (QR->QV)
       do k = KS, KE
          vent_factor = 1.6_RP + 124.9_RP * ( dens(k)*qr(k) )**0.2046_RP

          dq_evap(k) = ( 1.0_RP-min(Sliq(k),1.0_RP) ) / dens(k) * vent_factor  &
                     * ( dens(k)*qr(k) )**0.525_RP / ( 5.4E5_RP + 2.55E8_RP / ( pres(k)*qsatl(k) ) )
       enddo

       ! limiter
       do k = KS, KE
          dqc = (-dq_auto(k)-dq_accr(k)            )

          QTRC_t(k,i,j,I_QC) = QTRC_t(k,i,j,I_QC) + max( dqc, -qc(k)*rdt )
       enddo

       do k = KS, KE
          dqr = ( dq_auto(k)+dq_accr(k)-dq_evap(k) )

          QTRC_t(k,i,j,I_QR) = QTRC_t(k,i,j,I_QR) + max( dqr, -qr(k)*rdt )
       enddo

       do k = KS, KE
          dqv = - ( QTRC_t(k,i,j,I_QC) &
                  + QTRC_t(k,i,j,I_QR) )

          QTRC_t(k,i,j,I_QV) = QTRC_t(k,i,j,I_QV) + max( dqv, -qv(k)*rdt )
       enddo

       do k = KS, KE
          RHOE_t(k,i,j) = RHOE_t(k,i,j) - dens(k) * ( LHV00 * QTRC_t(k,i,j,I_QV) )
       enddo
    enddo
    enddo

    ! mass & energy update
    !$omp parallel do private(i,j, k) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
       QTRC0(KS:KE,i,j,I_QV) = QTRC0(KS:KE,i,j,I_QV) + QTRC_t(KS:KE,i,j,I_QV) * dt
       QTRC0(KS:KE,i,j,I_QC) = QTRC0(KS:KE,i,j,I_QC) + QTRC_t(KS:KE,i,j,I_QC) * dt
       QTRC0(KS:KE,i,j,I_QR) = QTRC0(KS:KE,i,j,I_QR) + QTRC_t(KS:KE,i,j,I_QR) * dt

       RHOE0(KS:KE,i,j) = RHOE0(KS:KE,i,j) + RHOE_t(KS:KE,i,j) * dt
    enddo
    enddo

    call PROF_rapend  ('MP_kessler')

    return
  end subroutine MP_kessler

  !-----------------------------------------------------------------------------
  !> Kessler-type warm rain microphysics (terminal velocity)
  !-----------------------------------------------------------------------------
  subroutine MP_kessler_vterm( &
       vterm, &
       DENS0, &
       QTRC0  )
    implicit none

    real(RP), intent(inout) :: vterm(KA,IA,JA,QA)
    real(RP), intent(in)    :: DENS0(KA,IA,JA)
    real(RP), intent(in)    :: QTRC0(KA,IA,JA,QA)

    real(RP) :: zerosw

    integer :: k, i, j
    !---------------------------------------------------------------------------

    ! only update QR
    !$omp parallel do private(i,j,k,zerosw) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       zerosw = 0.5_RP - sign(0.5_RP, QTRC0(k,i,j,I_QR) - 1.E-12_RP )
       vterm(k,i,j,I_QR) = - 36.34_RP * ( DENS0(k,i,j) * ( QTRC0(k,i,j,I_QR) + zerosw ) )**0.1364_RP &
                         * factor_vterm(k) * ( 1.0_RP - zerosw )
    enddo
    enddo
    enddo

    return
  end subroutine MP_kessler_vterm

  !-----------------------------------------------------------------------------
  !> Calculate Cloud Fraction
  subroutine ATMOS_PHY_MP_kessler_CloudFraction( &
       cldfrac, &
       QTRC     )
    use scale_tracer, only: &
       QAD => QA
    use scale_const, only: &
       EPS => CONST_EPS
    implicit none

    real(RP), intent(out) :: cldfrac(KA,IA,JA)
    real(RP), intent(in)  :: QTRC   (KA,IA,JA,QAD)

    real(RP) :: qhydro
    integer  :: k, i, j, iq
    !---------------------------------------------------------------------------

    do j  = JS, JE
    do i  = IS, IE
    do k  = KS, KE
       qhydro = 0.0_RP
       do iq = 1, MP_QA
          qhydro = qhydro + QTRC(k,i,j,I_MP2ALL(iq))
       enddo
       cldfrac(k,i,j) = 0.5_RP + sign(0.5_RP,qhydro-EPS)
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_MP_kessler_CloudFraction

  !-----------------------------------------------------------------------------
  !> Calculate Effective Radius
  subroutine ATMOS_PHY_MP_kessler_EffectiveRadius( &
       Re,    &
       QTRC0, &
       DENS0  )
    use scale_tracer, only: &
       QAD => QA, &
       MP_QAD => MP_QA
    implicit none

    real(RP), intent(out) :: Re   (KA,IA,JA,MP_QAD) ! effective radius
    real(RP), intent(in)  :: QTRC0(KA,IA,JA,QAD)    ! tracer mass concentration [kg/kg]
    real(RP), intent(in)  :: DENS0(KA,IA,JA)       ! density                   [kg/m3]
    !---------------------------------------------------------------------------

    Re(:,:,:,I_mp_QC) =   8.E-6_RP
    Re(:,:,:,I_mp_QR) = 100.E-6_RP

    return
  end subroutine ATMOS_PHY_MP_kessler_EffectiveRadius
  !-----------------------------------------------------------------------------
  !> Calculate mixing ratio of each category
  subroutine ATMOS_PHY_MP_kessler_Mixingratio( &
       Qe,    &
       QTRC0  )
    use scale_tracer, only: &
       QAD => QA, &
       MP_QAD => MP_QA
    use scale_const, only: &
       EPS => CONST_EPS
    implicit none

    real(RP), intent(out) :: Qe   (KA,IA,JA,MP_QAD) ! mixing ratio of each cateory [kg/kg]
    real(RP), intent(in)  :: QTRC0(KA,IA,JA,QAD)    ! tracer mass concentration [kg/kg]
    !---------------------------------------------------------------------------

    Qe(:,:,:,I_mp_QC) = QTRC0(:,:,:,I_QC)
    Qe(:,:,:,I_mp_QR) = QTRC0(:,:,:,I_QR)

    return
  end subroutine ATMOS_PHY_MP_kessler_Mixingratio
  !-----------------------------------------------------------------------------
end module scale_atmos_phy_mp_kessler
