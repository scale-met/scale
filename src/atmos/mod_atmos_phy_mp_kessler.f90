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
module mod_atmos_phy_mp
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_stdio, only: &
     IO_FID_LOG,  &
     IO_L
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_MP_setup
  public :: ATMOS_PHY_MP
  public :: MP_kessler
  public :: MP_kessler_vterm

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
  real(RP), public, save :: vterm(KA,IA,JA,QA) ! terminal velocity of each tracer [m/s]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical, private, save  :: MP_doreport_tendency = .false. ! report tendency of each process?
  logical, private, save  :: MP_donegative_fixer  = .true.  ! apply negative fixer?

  real(RP), private, save :: factor_vterm(KA) ! collection factor for terminal velocity of QR

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup Cloud Microphysics
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    use mod_comm, only: &
       COMM_horizontal_mean
    use mod_atmos_vars, only: &
       ATMOS_TYPE_PHY_MP, &
       DENS
    implicit none

    NAMELIST / PARAM_ATMOS_PHY_MP / &
       MP_doreport_tendency, &
       MP_donegative_fixer

    real(RP) :: rho_prof(KA) ! averaged profile of rho

    integer :: ierr
    integer :: k
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Cloud Microphisics]/Categ[ATMOS]'
    if( IO_L ) write(IO_FID_LOG,*) '*** KESSLER-type parametarization'

    if ( ATMOS_TYPE_PHY_MP /= 'KESSLER' ) then
       if ( IO_L ) write(IO_FID_LOG,*) 'xxx ATMOS_TYPE_PHY_MP is not KESSLER. Check!'
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

    ! Calculate collection factor for terminal velocity of QR
    call COMM_horizontal_mean( rho_prof(:), DENS(:,:,:) )
    rho_prof(:) = rho_prof(:) * 1.E-3_RP ! [kg/m3]->[g/cc]

    do k = KS, KE
       factor_vterm(k) = sqrt( rho_prof(KS)/rho_prof(k) )
    enddo

    vterm(:,:,:,:) = 0.0_RP

    return
  end subroutine ATMOS_PHY_MP_setup

  !-----------------------------------------------------------------------------
  !> Cloud Microphysics
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP
    use mod_atmos_vars, only: &
       ATMOS_vars_fillhalo, &
       ATMOS_vars_total,    &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       RHOT, &
       QTRC
    use mod_time, only: &
       dt => TIME_DTSEC_ATMOS_PHY_MP
    use mod_history, only: &
       HIST_in
    use mod_atmos_phy_mp_common, only: &
       MP_negative_fixer        => ATMOS_PHY_MP_negative_fixer,       &
       MP_precipitation         => ATMOS_PHY_MP_precipitation,        &
       MP_saturation_adjustment => ATMOS_PHY_MP_saturation_adjustment
    use mod_atmos_thermodyn, only: &
       THERMODYN_rhoe        => ATMOS_THERMODYN_rhoe,       &
       THERMODYN_rhot        => ATMOS_THERMODYN_rhot,       &
       THERMODYN_temp_pres_E => ATMOS_THERMODYN_temp_pres_E
    implicit none

    real(RP) :: RHOE_t(KA,IA,JA)
    real(RP) :: QTRC_t(KA,IA,JA,QA)
    real(RP) :: RHOE  (KA,IA,JA)
    real(RP) :: temp  (KA,IA,JA)
    real(RP) :: pres  (KA,IA,JA)

    real(RP) :: flux_tot (KA,IA,JA)
    real(RP) :: flux_rain(KA,IA,JA)
    real(RP) :: flux_snow(KA,IA,JA)
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Microphysics(kessler)'

    if ( MP_donegative_fixer ) then
       call MP_negative_fixer( DENS(:,:,:),  & ! [INOUT]
                               RHOT(:,:,:),  & ! [INOUT]
                               QTRC(:,:,:,:) ) ! [INOUT]
    endif

    call THERMODYN_rhoe( RHOE(:,:,:),  & ! [OUT]
                         RHOT(:,:,:),  & ! [IN]
                         QTRC(:,:,:,:) ) ! [IN]

    call MP_kessler( RHOE_t(:,:,:),   & ! [OUT]
                     QTRC_t(:,:,:,:), & ! [OUT]
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
                           temp (:,:,:)      ) ! [IN]

    call THERMODYN_rhot( RHOT(:,:,:),  & ! [OUT]
                         RHOE(:,:,:),  & ! [IN]
                         QTRC(:,:,:,:) ) ! [IN]

    call MP_saturation_adjustment( RHOT(:,:,:),   & ! [INOUT]
                                   QTRC(:,:,:,:), & ! [INOUT]
                                   DENS(:,:,:)    ) ! [IN]

    if ( MP_donegative_fixer ) then
       call MP_negative_fixer( DENS(:,:,:),  & ! [INOUT]
                               RHOT(:,:,:),  & ! [INOUT]
                               QTRC(:,:,:,:) ) ! [INOUT]
    endif

    ! fill halo
    call ATMOS_vars_fillhalo

    ! log report total (optional)
    call ATMOS_vars_total

    if ( MP_doreport_tendency ) then
       call HIST_in( QTRC_t(:,:,:,I_QV), 'QV_t_mp', 'tendency of QV in mp', 'kg/kg/s', dt )
       call HIST_in( QTRC_t(:,:,:,I_QC), 'QC_t_mp', 'tendency of QC in mp', 'kg/kg/s', dt )
       call HIST_in( QTRC_t(:,:,:,I_QR), 'QR_t_mp', 'tendency of QR in mp', 'kg/kg/s', dt )

       call HIST_in( RHOE_t(:,:,:), 'RHOE_t_mp', 'tendency of rhoe in mp', 'J/m3/s', dt )

       call HIST_in( vterm(:,:,:,I_QR), 'Vterm_QR', 'terminal velocity of QR', 'm/s', dt )
    endif

    flux_tot(:,:,:) = flux_rain(:,:,:) + flux_snow(:,:,:)
    call HIST_in( flux_rain(KS-1,:,:), 'RAIN', 'surface rain rate', 'kg/m2/s', dt)
    call HIST_in( flux_snow(KS-1,:,:), 'SNOW', 'surface snow rate', 'kg/m2/s', dt)
    call HIST_in( flux_tot (KS-1,:,:), 'PREC', 'surface precipitation rate', 'kg/m2/s', dt)

    return
  end subroutine ATMOS_PHY_MP

  !-----------------------------------------------------------------------------
  !> Kessler-type warm rain microphysics
  !-----------------------------------------------------------------------------
  subroutine MP_kessler( &
       RHOE_t, &
       QTRC_t, &
       RHOE0,  &
       QTRC0,  &
       DENS0   )
    use mod_const, only: &
       LHV00 => CONST_LH00
    use mod_time, only: &
       dt => TIME_DTSEC_ATMOS_PHY_MP
    use mod_atmos_thermodyn, only: &
       THERMODYN_temp_pres_E => ATMOS_THERMODYN_temp_pres_E
    use mod_atmos_saturation, only: &
       SATURATION_dens2qsat_liq => ATMOS_SATURATION_dens2qsat_liq
    implicit none

    real(RP), intent(out)   :: RHOE_t(KA,IA,JA)    ! tendency rhoe             [J/m3/s]
    real(RP), intent(out)   :: QTRC_t(KA,IA,JA,QA) ! tendency tracer           [kg/kg/s]
    real(RP), intent(inout) :: RHOE0 (KA,IA,JA)    ! density * internal energy [J/m3]
    real(RP), intent(inout) :: QTRC0 (KA,IA,JA,QA) ! mass concentration        [kg/kg]
    real(RP), intent(in)    :: DENS0 (KA,IA,JA)    ! density                   [kg/m3]

    ! working
    real(RP) :: QSATL(KA,IA,JA) ! saturated water vapor [kg/kg]
    real(RP) :: TEMP0(KA,IA,JA) ! temperature           [K]
    real(RP) :: PRES0(KA,IA,JA) ! pressure              [Pa]

    real(RP) :: dens
    real(RP) :: rhoe
    real(RP) :: temp
    real(RP) :: pres
    real(RP) :: q(QA)

    ! tendency
    real(RP) :: dq_evap ! tendency q (evaporation)
    real(RP) :: dq_auto ! tendency q (autoconversion)
    real(RP) :: dq_accr ! tendency q (accretion)
    real(RP) :: dqv, dqc, dqr
    real(RP) :: vent_factor, Sliq

    real(RP) :: rdt

    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    rdt = 1.D0 / dt

    call THERMODYN_temp_pres_E( TEMP0(:,:,:),  & ! [OUT]
                                PRES0(:,:,:),  & ! [OUT]
                                DENS0(:,:,:),  & ! [IN]
                                RHOE0(:,:,:),  & ! [IN]
                                QTRC0(:,:,:,:) ) ! [IN]

    call SATURATION_dens2qsat_liq( QSATL(:,:,:), & ! [OUT]
                                   TEMP0(:,:,:), & ! [IN]
                                   DENS0(:,:,:)  ) ! [IN]

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       ! store to work
       dens = DENS0(k,i,j)
       rhoe = RHOE0(k,i,j)
       temp = TEMP0(k,i,j)
       pres = PRES0(k,i,j)
       do iq = I_QV, I_QR
          q(iq) = QTRC0(k,i,j,iq)
       enddo

       Sliq = q(I_QV) / QSATL(k,i,j)

       ! Auto-conversion (QC->QR)
       dq_auto = 1.E-3_RP * max( q(I_QC)-1.E-3_RP, 0.0_RP )

       ! Accretion (QC->QR)
       dq_accr = 2.2_RP * q(I_QC) * q(I_QR)**0.875_RP

       ! Evaporation (QR->QV)
       vent_factor = 1.6_RP + 124.9_RP * ( dens*q(I_QR) )**0.2046_RP

       dq_evap = ( 1.0_RP-min(Sliq,1.0_RP) ) / dens * vent_factor  &
               * ( dens * q(I_QR) )**0.525_RP / ( 5.4E5_RP + 2.55E8_RP / ( pres*QSATL(k,i,j) ) )

       ! limiter
       dqc = ( -dq_auto-dq_accr         )
       dqr = (  dq_auto+dq_accr-dq_evap )

       ! tendency
       QTRC_t(k,i,j,I_QC) = max( dqc, -q(I_QC)*rdt )
       QTRC_t(k,i,j,I_QR) = max( dqr, -q(I_QR)*rdt )

       dqv = - ( QTRC_t(k,i,j,I_QC) &
               + QTRC_t(k,i,j,I_QR) )

       QTRC_t(k,i,j,I_QV) = max( dqv, -q(I_QV)*rdt )

       RHOE_t(k,i,j) = -dens * ( LHV00 * QTRC_t(k,i,j,I_QV) )
    enddo
    enddo
    enddo

    ! mass & energy update
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       QTRC0(k,i,j,I_QV) = QTRC0(k,i,j,I_QV) + QTRC_t(k,i,j,I_QV) * dt
       QTRC0(k,i,j,I_QC) = QTRC0(k,i,j,I_QC) + QTRC_t(k,i,j,I_QC) * dt
       QTRC0(k,i,j,I_QR) = QTRC0(k,i,j,I_QR) + QTRC_t(k,i,j,I_QR) * dt

       RHOE0(k,i,j) = RHOE0(k,i,j) + RHOE_t(k,i,j) * dt
    enddo
    enddo
    enddo

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

end module mod_atmos_phy_mp 
